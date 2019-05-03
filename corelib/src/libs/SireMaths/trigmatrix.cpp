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

#include "trigmatrix.h"
#include "nmatrix.h"
#include "matrix.h"
#include "nvector.h"
#include "vector.h"

#include "sire_blas.h"
#include "sire_lapack.h"

#include "SireBase/array2d.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<TrigMatrix> r_trigmatrix(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const TrigMatrix &matrix)
{
    writeHeader(ds, r_trigmatrix, 1);
    
    SharedDataStream sds(ds);
    
    sds << matrix.array << matrix.nrows;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, TrigMatrix &matrix)
{
    VersionID v = readHeader(ds, r_trigmatrix);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> matrix.array >> matrix.nrows;
    }
    else
        throw version_error(v, "1", r_trigmatrix, CODELOC);
        
    return ds;
}

/** Null constructor */
TrigMatrix::TrigMatrix() : nrows(0)
{}

/** Construct a matrix with dimension 'dimension*dimension'
    The values in the matrix are not initialised */
TrigMatrix::TrigMatrix(int dimension) : nrows(dimension)
{
    if (nrows <= 0)
    {
        nrows = 0;
    }
    else
    {
        array = QVector<double>( (nrows*nrows + nrows)/2 );
        array.squeeze();
    }
}

/** Construct a matrix with dimension 'dimension*dimension'
    The values in the matrix are initialised to 'initial_value' */
TrigMatrix::TrigMatrix(int dimension, double initial_value) : nrows(dimension)
{
    if (nrows <= 0)
    {
        nrows = 0;
    }
    else
    {
        array = QVector<double>( (nrows*nrows + nrows)/2, initial_value );
        array.squeeze();
    }
}

/** Construct from the passed Matrix - if 'take_upper' is true,
    then the upper right diagonal is copied - otherwise the 
    lower left diagonal is taken */
TrigMatrix::TrigMatrix(const NMatrix &matrix, bool take_upper) : nrows(0)
{
    if (matrix.nRows() > 0)
    {
        matrix.assertSquare();
        
        nrows = matrix.nRows();

        array = QVector<double>( (nrows*nrows + nrows)/2 );
        array.squeeze();
    
        double *data = array.data();
        const double *mdata = matrix.constData();
        
        if (take_upper)
        {
            for (int i=0; i<nrows; ++i)
            {
                for (int j=i; j<nrows; ++j)
                {
                    data[ offset(i,j) ] = mdata[ matrix.offset(i,j) ];
                }
            }
        }
        else
        {
            for (int i=0; i<nrows; ++i)
            {
                for (int j=i; j<nrows; ++j)
                {
                    data[ offset(j,i) ] = mdata[ matrix.offset(j,i) ];
                }
            }
        }
    }
}

/** Copy constructor */
TrigMatrix::TrigMatrix(const TrigMatrix &other)
        : array(other.array), nrows(other.nrows)
{}

/** Destructor */
TrigMatrix::~TrigMatrix()
{}

const char* TrigMatrix::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TrigMatrix>() );
}

const char* TrigMatrix::what() const
{
    return TrigMatrix::typeName();
}

/** Copy assignment operator */
TrigMatrix& TrigMatrix::operator=(const TrigMatrix &other)
{
    if (this != &other)
    {
        array = other.array;
        nrows = other.nrows;
    }

    return *this;
}

/** Comparison operator */
bool TrigMatrix::operator==(const TrigMatrix &other) const
{
    return this == &other or
           (nrows == other.nrows and array == other.array);
}

/** Comparison operator */
bool TrigMatrix::operator!=(const TrigMatrix &other) const
{
    return not this->operator==(other);
}

/** Assert that the index [i,j] is valid for this matrix

    \throw SireError::invalid_index
*/
void TrigMatrix::assertValidIndex(int i, int j) const
{
    if ( i < 0 or i >= nrows or j < 0 or j >= nrows )
    {
        throw SireError::invalid_index( QObject::tr(
            "The index [%1,%2] is invalid for the matrix with dimension [%1,%2].")
                .arg(i).arg(j).arg(nrows).arg(nrows), CODELOC );
    }
}

/** Assert that this is a square matrix - it definitely is! */
void TrigMatrix::assertSquare() const
{}

/** Return a modifiable reference to the value at [i,j]

    \throw SireError::invalid_index
*/
double& TrigMatrix::operator()(int i, int j)
{
    this->assertValidIndex(i,j);
    
    return array.data()[ this->offset(i,j) ];
}

/** Return a reference to the value at [i,j]

    \throw SireError::invalid_index
*/
const double& TrigMatrix::operator()(int i, int j) const
{
    this->assertValidIndex(i,j);
    
    return array.data()[ this->offset(i,j) ];
}

/** Assert that this matrix has 'nrows' rows

    \throw SireError::incompatible_error
*/
void TrigMatrix::assertNRows(int nr) const
{
    if (nr != nrows)
        throw SireError::incompatible_error( QObject::tr(
                "The number of rows in this matrix (dimension [%1,%2]) is not "
                "equal to %3.")
                    .arg(nrows).arg(nrows).arg(nr), CODELOC );
}

/** Assert that this matrix has 'ncolumns' columns

    \throw SireError::incompatible_error
*/
void TrigMatrix::assertNColumns(int nc) const
{
    if (nc != nrows)
        throw SireError::incompatible_error( QObject::tr(
                "The number of columns in this matrix (dimension [%1,%2]) is not "
                "equal to %3.")
                    .arg(nrows).arg(nrows).arg(nc), CODELOC );
}

/** Return the transpose of this matrix. This is fast, as it doesn't 
    actually have to do anything! */
TrigMatrix TrigMatrix::transpose() const
{
    return *this;
}

/** Fully transpose the data of this matrix - again, this does nothing */
TrigMatrix TrigMatrix::fullTranspose() const
{
    return *this;
}

/** Matrix addition 

    \throw SireError::incompatible_error
*/
TrigMatrix& TrigMatrix::operator+=(const TrigMatrix &other)
{
    assertNRows(other.nRows());
    
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
TrigMatrix& TrigMatrix::operator-=(const TrigMatrix &other)
{
    assertNRows(other.nRows());
    
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
TrigMatrix& TrigMatrix::operator*=(double scale)
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
TrigMatrix& TrigMatrix::operator/=(double scale)
{
    if (scale == 0)
        throw SireMaths::domain_error( QObject::tr(
            "This code does not support dividing a matrix by zero!"), CODELOC );
            
    return this->operator*=( 1/scale );
}

/** Matrix multiplication - this uses dgemm under the hood
    for speed 
    
    \throw SireError::incompatible_error
*/  
TrigMatrix& TrigMatrix::operator*=(const TrigMatrix &other)
{
    this->operator=( TrigMatrix( NMatrix(*this) * NMatrix(other) ) );
    
    return *this;
}

/** Return the inverse of this matrix
    
    This uses LAPACK under the hood, for speed
    
    \throw SireError::incompatible_error
    \throw SireMaths::domain_error
*/
TrigMatrix TrigMatrix::inverse() const
{
    return TrigMatrix( NMatrix(*this).inverse() );
}

/** Matrix division - this multiplies this matrix with the inverse of 'other' 

    \throw SireMaths::domain_error
    \throw SireError::incompatible_error
*/
TrigMatrix& TrigMatrix::operator/=(const TrigMatrix &other)
{
    this->operator=( TrigMatrix( NMatrix(*this) / NMatrix(other) ) );
    return *this;
}

/** Return the negative of this matrix */
TrigMatrix TrigMatrix::operator-() const
{
    TrigMatrix ret(*this);

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
TrigMatrix TrigMatrix::operator+(const TrigMatrix &other) const
{
    TrigMatrix ret(*this);
    ret += other;
    
    return ret;
}

/** Matrix subtraction 

    \throw SireError::incompatible_error
*/
TrigMatrix TrigMatrix::operator-(const TrigMatrix &other) const
{
    TrigMatrix ret(*this);
    ret -= other;
    
    return ret;
}

/** Matrix multiplication - this uses dgemm under the hood
    for speed 
    
    \throw SireError::incompatible_error
*/  
TrigMatrix TrigMatrix::operator*(const TrigMatrix &other) const
{
    TrigMatrix ret(*this);
    
    ret *= other;
    
    return ret;
}

/** Matrix division - this multiplies this matrix with the inverse of 'other' 

    \throw SireMaths::domain_error
    \throw SireError::incompatible_error
*/
TrigMatrix TrigMatrix::operator/(const TrigMatrix &other) const
{
    TrigMatrix ret(*this);
    
    ret /= other;
    
    return ret;
}

/** Multiply all elements of this matrix by 'scale' */
TrigMatrix TrigMatrix::operator*(double scale) const
{
    TrigMatrix ret(*this);
    
    ret *= scale;
    
    return ret;
}

/** Divide all elements of this matrix by 'scale' */
TrigMatrix TrigMatrix::operator/(double scale) const
{
    TrigMatrix ret(*this);
    
    ret /= scale;
    
    return ret;
}

/** Perform matrix-vector multiplication - the number of 
    rows of the vector must equal to the number of rows
    of this matrix

    This uses dgemv under the hood for speed

    \throw SireError::incompatible_error
*/
NVector TrigMatrix::operator*(const NVector &vector) const
{
    return NMatrix(*this) * vector;
}

/** Perform matrix-vector multiplication - the number of 
    rows of the vector must equal to the number of rows
    of this matrix

    This uses dgemv under the hood for speed

    \throw SireError::incompatible_error
*/
NVector TrigMatrix::operator*(const Vector &vector) const
{
    return NMatrix(*this) * vector;
}

/** Return the number of rows in this matrix */
int TrigMatrix::nRows() const
{
    return nrows;
}

/** Return the number of columns in this matrix */
int TrigMatrix::nColumns() const
{
    return nrows;
}

/** Redimension this matrix to have 'dimension' rows and 'dimension' 
    columns. The contents of this matrix are undefined after
    this redimension. This function will only reallocate
    memory if there is not enough memory allocated to store
    the new matrix. Use this function if you want to use
    the same piece of memory over and over again for lots
    of different size matricies - just create a matrix with
    the maximum dimension, then call this 'redimension' function
    whenever you want to change. It is very fast, as it just
    updates the internal record of the size of the matrix */
void TrigMatrix::redimension(int dimension)
{
    const int sz = (dimension*dimension + dimension) / 2;
    
    if (sz <= 0)
    {
        nrows = 0;
    }
    else
    {
        if (sz > array.count())
            array.resize(sz);
            
        nrows = dimension;
    }
}

/** Assert that there is an ith row! 

    \throw SireError::invalid_index
*/
void TrigMatrix::assertValidRow(int i) const
{
    if (i < 0 or i >= nrows)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2] does not have a row with index %3.")
                    .arg(nrows).arg(nrows).arg(i), CODELOC );
}

/** Assert that there is an jth column! 

    \throw SireError::invalid_index
*/
void TrigMatrix::assertValidColumn(int j) const
{
    if (j < 0 or j >= nrows)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2] does not have a column with index %3.")
                    .arg(nrows).arg(nrows).arg(j), CODELOC );
}

/** Return a vector containing the contents of the ith row 

    \throw SireError::invalid_index
*/
NVector TrigMatrix::row(int i) const
{
    this->assertValidRow(i);

    NVector v(nrows);

    const double *d = array.constData();
    double *vd = v.data();
    
    for (int j=0; j<nrows; ++j)
    {
        vd[j] = d[ offset(i,j) ];
    }

    return v;
}

/** Return a vector containing the contents of the jth column

    \throw SireError::invalid_index
*/
NVector TrigMatrix::column(int j) const
{
    this->assertValidColumn(j);

    NVector v(nrows);
    const double *d = array.constData();
    double *vd = v.data();
    
    for (int i=0; i<nrows; ++i)
    {
        vd[i] = d[ offset(i,j) ];
    }

    return v;
}

/** Set the value of [i,j] (and [j,i]) to 'value'

    \throw SireError::invalid_index
*/
void TrigMatrix::set(int i, int j, double value)
{
    array.data()[checkedOffset(i,j)] = value;
}

/** Set the values of all data in the row 'i' to 'value'

    \throw SireError::invalid_index
*/
void TrigMatrix::setRow(int i, double value)
{
    this->assertValidRow(i);
    
    double *d = array.data();
    
    for (int j=0; j<nrows; ++j)
    {
        d[offset(i,j)] = value;
    }
}

/** Copy the vector 'row' to row 'i'

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void TrigMatrix::setRow(int i, const NVector &row)
{
    this->assertValidRow(i);
    this->assertNColumns(row.count());
    
    double *d = array.data();
    const double *v = row.constData();
    
    for (int j=0; j<nrows; ++j)
    {
        d[offset(i,j)] = v[j];
    }
}

/** Set the values of all data in the column 'j' to 'value'

    \throw SireError::invalid_index
*/
void TrigMatrix::setColumn(int j, double value)
{
    this->assertValidColumn(j);
    
    double *d = array.data();
    
    for (int i=0; i<nrows; ++i)
    {
        d[offset(i,j)] = value;
    }
}

/** Copy the vector 'column' to column 'j'

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void TrigMatrix::setColumn(int j, const NVector &column)
{
    this->assertValidColumn(j);
    this->assertNRows(column.count());
    
    double *d = array.data();
    const double *v = column.constData();
    
    for (int i=0; i<nrows; ++i)
    {
        d[offset(i,j)] = v[i];
    }
}

/** Set all entries in the matrix to the value 'value' */
void TrigMatrix::setAll(double value)
{
    double *d = array.data();
    int sz = array.count();
    
    for (int i=0; i<sz; ++i)
    {
        d[i] = value;
    }
}

/** Return a raw pointer to the data of this matrix. The data
    is stored in a packed format - use 'offset' or 'checkedOffset'
    to get the index of the element 'i,j' */
double* TrigMatrix::data()
{
    return array.data();
}

/** Return a raw pointer to the data of this matrix. The data
    is stored in a packed format - use 'offset' or 'checkedOffset'
    to get the index of the element 'i,j' */
const double* TrigMatrix::data() const
{
    return array.constData();
}

/** Return a raw pointer to the data of this matrix. The data
    is stored in a packed format - use 'offset' or 'checkedOffset'
    to get the index of the element 'i,j' */
const double* TrigMatrix::constData() const
{
    return array.constData();
}

/** Return the QVector containing the memory of this Matrix */
QVector<double> TrigMatrix::memory() const
{
    return array;
}

/** Calculate the offset in the 1D array of the value
    at index [i,j]
    
    \throw SireError::invalid_index
*/
int TrigMatrix::checkedOffset(int i, int j) const
{
    this->assertValidIndex(i,j);
    return this->offset(i,j);
}

/** Return a string representation of this matrix */
QString TrigMatrix::toString() const
{
    if (nrows == 0)
        return "( )";
        
    else if (nrows == 1)
    {
        const double *d = array.constData();
    
        QStringList row;
        for (qint32 j=0; j<nrows; ++j)
        {
            row.append( QString("%1").arg(d[j], 8) );
        }
        
        return QString("( %1 )").arg( row.join(", ") );
    }

    QStringList rows;
    
    const double *d = array.constData();
    
    for (qint32 i=0; i<nrows; ++i)
    {
        QStringList row;
        
        for (qint32 j=0; j<nrows; ++j)
        {
            row.append( QString("%1").arg(d[offset(i,j)], 8) );
        }
        
        if (i == 0)
            rows.append( QString("/ %1 \\").arg( row.join(", ") ) );
        else if (i == this->nRows() - 1)
            rows.append( QString("\\ %1 /").arg( row.join(", ") ) );
        else
            rows.append( QString("| %1 |").arg( row.join(", ") ) );
    }
    
    return rows.join("\n");
}

/** Return the determinant of this matrix

    This uses LAPACK under the hood, for speed

    \throw SireError::incompatible_error
*/
double TrigMatrix::determinant() const
{
    return NMatrix(*this).determinant();
}

/** Return the trace of this matrix - this is only valid for a square matrix

    \throw SireError::incompatible_error
*/
double TrigMatrix::trace() const
{
    const double *d = array.constData();
    double sum = 0;
    
    for (int i=0; i<nrows; ++i)
    {
        sum += d[ offset(i,i) ];
    }
    
    return sum;
}
    
/** Return a vector containing the diagonal of this matrix - this is only
    valid for a square matrix
    
    \throw SireError::incompatible_error
*/
NVector TrigMatrix::diagonal() const
{
    const double *d = array.constData();

    NVector vector(nrows);
    double *v = vector.data();
    
    for (int i=0; i<nrows; ++i)
    {
        v[i] = d[ offset(i,i) ];
    }
    
    return vector;
}

/** Return whether or not this is a transposed matrix (data
    is stored in row-major order rather than column-major order) */
bool TrigMatrix::isTransposed() const
{
    return false;
}

/** Return the eigenvalues and eigenvectors of this matrix. This
    uses LAPACK under the hood for speed

    \throw SireError::incompatible_error
    \throw SireMaths::domain_error
*/
std::pair<NVector,NMatrix> TrigMatrix::diagonalise() const
{
    return NMatrix(*this).diagonalise();
}
