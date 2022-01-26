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

#include "nmatrix.h"
#include "trigmatrix.h"
#include "matrix.h"
#include "nvector.h"
#include "vector.h"

#include "sire_blas.h"
#include "sire_lapack.h"
#include "sire_linpack.h"

#include "SireBase/array2d.hpp"
#include "SireBase/trigarray2d.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<NMatrix> r_nmatrix(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NMatrix &matrix)
{
    writeHeader(ds, r_nmatrix, 1);

    SharedDataStream sds(ds);

    sds << matrix.array << matrix.nrows << matrix.ncolumns << matrix.is_transpose;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NMatrix &matrix)
{
    VersionID v = readHeader(ds, r_nmatrix);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> matrix.array >> matrix.nrows >> matrix.ncolumns >> matrix.is_transpose;
    }
    else
        throw version_error(v, "1", r_nmatrix, CODELOC);

    return ds;
}

/** Null constructor */
NMatrix::NMatrix() : nrows(0), ncolumns(0), is_transpose(false)
{}

/** Construct a matrix with 'nrows' rows and 'ncolumns' columns.
    The values in the matrix are not initialised */
NMatrix::NMatrix(int nr, int nc)
        : nrows(nr), ncolumns(nc), is_transpose(false)
{
    if (nr <= 0 or nc <= 0)
    {
        nrows = 0;
        ncolumns = 0;
    }
    else
    {
        array = QVector<double>(nrows*ncolumns);
        array.squeeze();
    }
}

/** Construct a matrix with 'nrows' rows and 'ncolumns' columns.
    The values in the matrix are initialised to 'initial_value' */
NMatrix::NMatrix(int nr, int nc, double initial_value)
        : nrows(nr), ncolumns(nc), is_transpose(false)
{
    if (nr <= 0 or nc <= 0)
    {
        nrows = 0;
        ncolumns = 0;
    }
    else
    {
        array = QVector<double>(nrows*ncolumns, initial_value);
        array.squeeze();
    }
}

/** Construct from the passed Matrix */
NMatrix::NMatrix(const Matrix &matrix)
        : nrows(3), ncolumns(3), is_transpose(false)
{
    array = QVector<double>(9);
    array.squeeze();

    double *data = array.data();

    memcpy(data, matrix.column0().constData(), 3*sizeof(double));
    memcpy(data+3, matrix.column1().constData(), 3*sizeof(double));
    memcpy(data+6, matrix.column2().constData(), 3*sizeof(double));
}

/** Construct from the passed Array */
NMatrix::NMatrix(const SireBase::Array2D<double> &array2d)
        : nrows(array2d.nRows()), ncolumns(array2d.nColumns()), is_transpose(true)
{
    int sz = array2d.nRows() * array2d.nColumns();

    if (sz > 0)
    {
        array = QVector<double>(sz);
        array.squeeze();

        memcpy(array.data(), array2d.constData(), sz*sizeof(double));
    }
}

/** Construct from the passed array */
NMatrix::NMatrix(const QVector< QVector<double> > &array2d)
        : nrows(array2d.count()), ncolumns(0), is_transpose(true)
{
    const QVector<double> *array2d_data = array2d.constData();

    for (int i=0; i<nrows; ++i)
    {
        ncolumns = qMax(ncolumns, array2d_data[i].count());
    }

    if (ncolumns > 0)
    {
        array = QVector<double>(nrows*ncolumns, 0);
        array.squeeze();

        double *d = array.data();

        for (int i=0; i<nrows; ++i)
        {
            const double *array_d = array2d_data[i].constData();
            const int sz = array2d_data[i].count();

            for (int j=0; j<sz; ++j)
            {
                d[offset(i,j)] = array_d[j];
            }
        }
    }
    else
    {
        nrows = 0;
    }
}

/** Construct from the passed vector - this is copied to a column
    matrix, unless 'transpose' is true, in which case this is
    copied as a row matrix */
NMatrix::NMatrix(const Vector &vector, bool transpose)
        : is_transpose(false)
{
    array = QVector<double>(3);
    array.squeeze();

    memcpy(array.data(), vector.constData(), 3*sizeof(double));

    if (transpose)
    {
        nrows = 1;
        ncolumns = 3;
    }
    else
    {
        nrows = 3;
        ncolumns = 1;
    }
}

/** Construct from the passed vector - this is copied to a column
    matrix, unless 'transpose' is true, in which case this is
    copied as a row matrix */
NMatrix::NMatrix(const NVector &vector, bool transpose)
        : nrows(0), ncolumns(0), is_transpose(false)
{
    if (vector.count() > 0)
    {
        array = QVector<double>(vector.count());
        array.squeeze();

        memcpy(array.data(), vector.constData(), vector.count()*sizeof(double));

        if (transpose)
        {
            nrows = 1;
            ncolumns = vector.count();
        }
        else
        {
            nrows = vector.count();
            ncolumns = 1;
        }
    }
}

/** Construct from the passed vector - this is copied to a column
    matrix, unless 'transpose' is true, in which case this is
    copied as a row matrix */
NMatrix::NMatrix(const QVector<double> &vector, bool transpose)
        : nrows(0), ncolumns(0), is_transpose(false)
{
    if (vector.count() > 0)
    {
        array = vector;
        array.squeeze();

        if (transpose)
        {
            nrows = 1;
            ncolumns = vector.count();
        }
        else
        {
            nrows = vector.count();
            ncolumns = 1;
        }
    }
}

/** Construct from the passed triangular matrix */
NMatrix::NMatrix(const TrigMatrix &matrix)
        : nrows(0), ncolumns(0), is_transpose(false)
{
    if (matrix.nRows() > 0)
    {
        nrows = matrix.nRows();
        ncolumns = nrows;

        array = QVector<double>(nrows * ncolumns);
        array.squeeze();

        double *d = array.data();
        const double *md = matrix.constData();

        for (int j=0; j<nrows; ++j)
        {
            for (int i=0; i<nrows; ++i)
            {
                *d = md[ matrix.offset(i,j) ];
                ++d;
            }
        }
    }
}

/** Copy constructor */
NMatrix::NMatrix(const NMatrix &other)
        : array(other.array), nrows(other.nrows),
          ncolumns(other.ncolumns), is_transpose(other.is_transpose)
{}

/** Destructor */
NMatrix::~NMatrix()
{}

/** Construct a matrix with dimension 'nrows' by 'ncolumns' that
    stores the data in column-major order */
NMatrix NMatrix::createColumnMajor(int nr, int nc)
{
    return NMatrix(nr, nc);
}

/** Construct a matrix with dimension 'nrows' by 'ncolumns' that
    stores the data in row-major order */
NMatrix NMatrix::createRowMajor(int nr, int nc)
{
    NMatrix ret(nr, nc);
    ret.is_transpose = true;
    return ret;
}

const char* NMatrix::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NMatrix>() );
}

const char* NMatrix::what() const
{
    return NMatrix::typeName();
}

/** Copy assignment operator */
NMatrix& NMatrix::operator=(const NMatrix &other)
{
    if (this != &other)
    {
        array = other.array;
        nrows = other.nrows;
        ncolumns = other.ncolumns;
        is_transpose = other.is_transpose;
    }

    return *this;
}

/** Comparison operator */
bool NMatrix::operator==(const NMatrix &other) const
{
    return this == &other or
           (nrows == other.nrows and ncolumns == other.ncolumns and
            is_transpose == other.is_transpose and
            array == other.array);
}

/** Comparison operator */
bool NMatrix::operator!=(const NMatrix &other) const
{
    return not this->operator==(other);
}

/** Assert that the index [i,j] is valid for this matrix

    \throw SireError::invalid_index
*/
void NMatrix::assertValidIndex(int i, int j) const
{
    if ( i < 0 or i >= nrows or j < 0 or j >= ncolumns )
    {
        throw SireError::invalid_index( QObject::tr(
            "The index [%1,%2] is invalid for the matrix with dimension [%1,%2].")
                .arg(i).arg(j).arg(nrows).arg(ncolumns), CODELOC );
    }
}

/** Assert that this is a square matrix

    \throw SireError::incompatible_error
*/
void NMatrix::assertSquare() const
{
    if (nrows == 0 or ncolumns == 0)
        throw SireError::incompatible_error( QObject::tr(
            "The null matrix of zero dimension is not a square matrix!"), CODELOC );

    else if (nrows != ncolumns)
        throw SireError::incompatible_error( QObject::tr(
            "The operation is only compatible with square matricies. This "
            "matrix, of dimension [%1,%2], is not a square matrix.")
                .arg(nrows).arg(ncolumns), CODELOC );
}

/** Return a modifiable reference to the value at [i,j]

    \throw SireError::invalid_index
*/
double& NMatrix::operator()(int i, int j)
{
    this->assertValidIndex(i,j);

    return array.data()[ this->offset(i,j) ];
}

/** Return a reference to the value at [i,j]

    \throw SireError::invalid_index
*/
const double& NMatrix::operator()(int i, int j) const
{
    this->assertValidIndex(i,j);

    return array.data()[ this->offset(i,j) ];
}

/** Assert that this matrix has 'nrows' rows

    \throw SireError::incompatible_error
*/
void NMatrix::assertNRows(int nr) const
{
    if (nr != nrows)
        throw SireError::incompatible_error( QObject::tr(
                "The number of rows in this matrix (dimension [%1,%2]) is not "
                "equal to %3.")
                    .arg(nrows).arg(ncolumns).arg(nr), CODELOC );
}

/** Assert that this matrix has 'ncolumns' columns

    \throw SireError::incompatible_error
*/
void NMatrix::assertNColumns(int nc) const
{
    if (nc != ncolumns)
        throw SireError::incompatible_error( QObject::tr(
                "The number of columns in this matrix (dimension [%1,%2]) is not "
                "equal to %3.")
                    .arg(nrows).arg(ncolumns).arg(nc), CODELOC );
}

/** Return the transpose of this matrix. This is fast, as this
    just toggles a flag to say whether or not the transpose is
    to be used. If you want to fully transpose the data (e.g.
    if you want to directly access the data) the call 'fullTranspose()' */
NMatrix NMatrix::transpose() const
{
    /*   / a b c \      / a d g j \
         | d e f |  =>  | b e h k |
         | g h i |      \ c f i l /
         \ j k l /

       in memory stored as column-major
       (a d g j b e h k c f i l)

        offset(i,j)   == i + (j*nrows)     [2,1] => 2+4 == 6 == 'h'
        offset(i,j)^T == (i*ncolumns) + j  [1,2] => 4+2 == 6 == 'h'
    */

    NMatrix ret(*this);
    ret.is_transpose = not is_transpose;
    ret.nrows = ncolumns;
    ret.ncolumns = nrows;

    return ret;
}

/** Fully transpose the data of this matrix */
NMatrix NMatrix::fullTranspose() const
{
    NMatrix ret(*this);
    ret.is_transpose = false;
    ret.ncolumns = nrows;
    ret.nrows = ncolumns;

    if (not is_transpose)
    {
        //the data needs to be transposed
        double *ret_data = ret.array.data();
        const double *this_data = array.constData();

        for (qint32 i=0; i<nrows; ++i)
        {
            for (qint32 j=0; j<ncolumns; ++j)
            {
                ret_data[ ret.offset(j,i) ] = this_data[ offset(i,j) ];
            }
        }
    }

    return ret;
}

/** Matrix addition

    \throw SireError::incompatible_error
*/
NMatrix& NMatrix::operator+=(const NMatrix &other)
{
    assertNRows(other.nRows());
    assertNColumns(other.nColumns());

    if (is_transpose == other.is_transpose)
    {
        double *data = array.data();
        const double *other_data = other.array.constData();
        const int sz = array.count();

        for (int i=0; i<sz; ++i)
        {
            data[i] += other_data[i];
        }
    }
    else
    {
        this->operator+=(other.fullTranspose());
    }

    return *this;
}

/** Matrix subtraction

    \throw SireError::incompatible_error
*/
NMatrix& NMatrix::operator-=(const NMatrix &other)
{
    assertNRows(other.nRows());
    assertNColumns(other.nColumns());

    if (is_transpose == other.is_transpose)
    {
        double *data = array.data();
        const double *other_data = other.array.constData();
        const int sz = array.count();

        for (int i=0; i<sz; ++i)
        {
            data[i] -= other_data[i];
        }
    }
    else
    {
        this->operator-=(other.fullTranspose());
    }

    return *this;
}

/** Multiply all elements of this matrix by 'scale' */
NMatrix& NMatrix::operator*=(double scale)
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
            data[i] *= scale;
        }
    }

    return *this;
}

/** Divide all elements of this matrix by 'scale' */
NMatrix& NMatrix::operator/=(double scale)
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
NMatrix& NMatrix::operator*=(const NMatrix &other)
{
    if (other.nRows() == 1 and other.nColumns() == 1)
        return this->operator*=(other(0,0));

    this->operator=( dgemm(*this, other) );

    return *this;
}

/** Return the inverse of this matrix

    This uses LAPACK under the hood, for speed

    \throw SireError::incompatible_error
    \throw SireMaths::domain_error
*/
NMatrix NMatrix::inverse() const
{
    this->assertSquare();
    std::pair< NMatrix,QVector<int> > factors_with_pivot = dgeco(*this);
    return dgedi_inverse(factors_with_pivot.first,
                         factors_with_pivot.second);
}

/** Matrix division - this multiplies this matrix with the inverse of 'other'

    \throw SireMaths::domain_error
    \throw SireError::incompatible_error
*/
NMatrix& NMatrix::operator/=(const NMatrix &other)
{
    return this->operator*=(other.inverse());
}

/** Return the negative of this matrix */
NMatrix NMatrix::operator-() const
{
    NMatrix ret(*this);

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
NMatrix NMatrix::operator+(const NMatrix &other) const
{
    NMatrix ret(*this);
    ret += other;

    return ret;
}

/** Matrix subtraction

    \throw SireError::incompatible_error
*/
NMatrix NMatrix::operator-(const NMatrix &other) const
{
    NMatrix ret(*this);
    ret -= other;

    return ret;
}

/** Matrix multiplication - this uses dgemm under the hood
    for speed

    \throw SireError::incompatible_error
*/
NMatrix NMatrix::operator*(const NMatrix &other) const
{
    NMatrix ret(*this);

    ret *= other;

    return ret;
}

/** Matrix division - this multiplies this matrix with the inverse of 'other'

    \throw SireMaths::domain_error
    \throw SireError::incompatible_error
*/
NMatrix NMatrix::operator/(const NMatrix &other) const
{
    NMatrix ret(*this);

    ret /= other;

    return ret;
}

/** Multiply all elements of this matrix by 'scale' */
NMatrix NMatrix::operator*(double scale) const
{
    NMatrix ret(*this);

    ret *= scale;

    return ret;
}

/** Divide all elements of this matrix by 'scale' */
NMatrix NMatrix::operator/(double scale) const
{
    NMatrix ret(*this);

    ret /= scale;

    return ret;
}

/** Perform matrix-vector multiplication - the number of
    rows of the vector must equal to the number of rows
    of this matrix

    This uses dgemv under the hood for speed

    \throw SireError::incompatible_error
*/
NVector NMatrix::operator*(const NVector &vector) const
{
    return dgemv(*this, vector);
}

/** Perform matrix-vector multiplication - the number of
    rows of the vector must equal to the number of rows
    of this matrix

    This uses dgemv under the hood for speed

    \throw SireError::incompatible_error
*/
NVector NMatrix::operator*(const Vector &vector) const
{
    if (this->nRows() == 3 and this->nColumns() == 3)
    {
        //use hand-written code
        const double *d = array.constData();

        if (is_transpose)
        {
            Matrix m( d[0], d[1], d[2],
                      d[3], d[4], d[5],
                      d[6], d[7], d[8] );   // Matrix uses row-major ordering

            return m*vector;
        }
        else
        {
            Matrix m( d[0], d[3], d[6],
                      d[1], d[4], d[7],
                      d[2], d[5], d[8] );   // Matrix uses row-major ordering

            return m*vector;
        }
    }
    else
    {
        this->assertNColumns(3);

        NVector v(nrows);

        const double *d = array.data();

        for (int i=0; i<nrows; ++i)
        {
            double sum = 0;

            for (int j=0; j<3; ++j)
            {
                sum += vector[j]*d[offset(i,j)];
            }
        }

        return v;
    }
}

/** Return the number of rows in this matrix */
int NMatrix::nRows() const
{
    return nrows;
}

/** Return the number of columns in this matrix */
int NMatrix::nColumns() const
{
    return ncolumns;
}

/** Redimension this matrix to have 'nrows' rows and 'ncolumns'
    columns. The contents of this matrix are undefined after
    this redimension. This function will only reallocate
    memory if there is not enough memory allocated to store
    the new matrix. Use this function if you want to use
    the same piece of memory over and over again for lots
    of different size matricies - just create a matrix with
    the maximum dimension, then call this 'redimension' function
    whenever you want to change. It is very fast, as it just
    updates the internal record of the size of the matrix */
void NMatrix::redimension(int nr, int nc)
{
    const int sz = nr * nc;

    if (sz <= 0)
    {
        nrows = 0;
        ncolumns = 0;
    }
    else
    {
        if (sz > array.count())
            array.resize(sz);

        nrows = nr;
        ncolumns = nc;

        is_transpose = false;
    }
}

/** Assert that there is an ith row!

    \throw SireError::invalid_index
*/
void NMatrix::assertValidRow(int i) const
{
    if (i < 0 or i >= nrows)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2] does not have a row with index %3.")
                    .arg(nrows).arg(ncolumns).arg(i), CODELOC );
}

/** Assert that there is an jth column!

    \throw SireError::invalid_index
*/
void NMatrix::assertValidColumn(int j) const
{
    if (j < 0 or j >= ncolumns)
        throw SireError::invalid_index( QObject::tr(
                "The matrix with dimension [%1,%2] does not have a column with index %3.")
                    .arg(nrows).arg(ncolumns).arg(j), CODELOC );
}

/** Return a vector containing the contents of the ith row

    \throw SireError::invalid_index
*/
NVector NMatrix::row(int i) const
{
    this->assertValidRow(i);

    NVector v(ncolumns);

    if (is_transpose)
    {
        //row-major storage
        memcpy( v.data(), array.constData()+i*ncolumns, ncolumns*sizeof(double) );
    }
    else
    {
        //column-major storage
        double *d = v.data();
        const double *row = array.constData();

        for (int j=0; j<ncolumns; ++j)
        {
            d[j] = row[ j*nrows + i ];
        }
    }

    return v;
}

/** Return a vector containing the contents of the ith row

    \throw SireError::invalid_index
*/
NVector NMatrix::column(int j) const
{
    this->assertValidColumn(j);

    NVector v(nrows);

    if (is_transpose)
    {
        //row-major storage
        double *d = v.data();
        const double *column = array.constData();

        for (int i=0; i<nrows; ++i)
        {
            d[i] = column[ i*ncolumns + j ];
        }
    }
    else
    {
        //column-major storage
        memcpy( v.data(), array.constData()+j*nrows, nrows*sizeof(double) );
    }

    return v;
}

/** Set the value of [i,j] to 'value'

    \throw SireError::invalid_index
*/
void NMatrix::set(int i, int j, double value)
{
    array.data()[checkedOffset(i,j)] = value;
}

/** Set the values of all data in the row 'i' to 'value'

    \throw SireError::invalid_index
*/
void NMatrix::setRow(int i, double value)
{
    this->assertValidRow(i);

    double *d = array.data();

    for (int j=0; j<ncolumns; ++j)
    {
        d[offset(i,j)] = value;
    }
}

/** Copy the vector 'row' to row 'i'

    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void NMatrix::setRow(int i, const NVector &row)
{
    this->assertValidRow(i);
    this->assertNColumns(row.count());

    double *d = array.data();
    const double *v = row.constData();

    if (is_transpose)
    {
        memcpy(d + offset(i,0), v, ncolumns*sizeof(double));
    }
    else
    {
        for (int j=0; j<ncolumns; ++j)
        {
            d[offset(i,j)] = v[j];
        }
    }
}

/** Set the values of all data in the column 'j' to 'value'

    \throw SireError::invalid_index
*/
void NMatrix::setColumn(int j, double value)
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
void NMatrix::setColumn(int j, const NVector &column)
{
    this->assertValidColumn(j);
    this->assertNRows(column.count());

    double *d = array.data();
    const double *v = column.constData();

    if (is_transpose)
    {
        for (int i=0; i<nrows; ++i)
        {
            d[offset(i,j)] = v[i];
        }
    }
    else
    {
        memcpy(d + offset(0,j), v, nrows*sizeof(double));
    }
}

/** Set all entries in the matrix to the value 'value' */
void NMatrix::setAll(double value)
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
double* NMatrix::data()
{
    return array.data();
}

/** Return a raw pointer to the data of this matrix. The data is
    stored in column-major order (same as Fortran - not same as C++ or C).
    To be safe, use the 'offset' function to get the offset of
    the value at [i,j] in this array */
const double* NMatrix::data() const
{
    return array.constData();
}

/** Return a raw pointer to the data of this matrix. The data is
    stored in column-major order (same as Fortran - not same as C++ or C).
    To be safe, use the 'offset' function to get the offset of
    the value at [i,j] in this array */
const double* NMatrix::constData() const
{
    return array.constData();
}

/** Return the QVector containing the memory of this Matrix */
QVector<double> NMatrix::memory() const
{
    return array;
}

/** Calculate the offset in the 1D array of the value
    at index [i,j]

    \throw SireError::invalid_index
*/
int NMatrix::checkedOffset(int i, int j) const
{
    this->assertValidIndex(i,j);
    return this->offset(i,j);
}

/** Return a string representation of this matrix */
QString NMatrix::toString() const
{
    if (nrows == 0)
        return "( )";

    else if (nrows == 1)
    {
        const double *d = array.constData();

        QStringList row;
        for (qint32 j=0; j<ncolumns; ++j)
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

        for (qint32 j=0; j<ncolumns; ++j)
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

/** Reflect the contents of the top half to the bottom
    half. If n == nRows(), then this sets

    matrix[n-i,j] = matrix[i,j]

    1 2 3      1 2 3
    4 5 6  =>  4 5 6
    7 8 9      1 2 3

    \throw SireError::incompatible_error
*/
void NMatrix::reflectTopToBottom()
{
    double *data = array.data();

    for (int i=0; i<nrows/2; ++i)
    {
        for (int j=0; j<ncolumns; ++j)
        {
            data[ offset(nrows-i-1,j) ] = data[ offset(i,j) ];
        }
    }
}

/** Reflect the contents of the bottom half to the top
    half. If n == nRows(), then this sets

    matrix[i,j] = matrix[n-i,j]

    1 2 3      7 8 9
    4 5 6  =>  4 5 6
    7 8 9      7 8 9

    \throw SireError::incompatible_error
*/
void NMatrix::reflectBottomToTop()
{
    double *data = array.data();

    for (int i=0; i<nrows/2; ++i)
    {
        for (int j=0; j<ncolumns; ++j)
        {
            data[ offset(i,j) ] = data[ offset(nrows-i-1,j) ];
        }
    }
}

/** Reflect the contents of the left half to the right
    half. If n == nColumns(), then this sets

    matrix[i,n-j] = matrix[i,j]

    1 2 3      1 2 1
    4 5 6  =>  4 5 4
    7 8 9      7 8 7

    \throw SireError::incompatible_error
*/
void NMatrix::reflectLeftToRight()
{
    double *data = array.data();

    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncolumns/2; ++j)
        {
            data[ offset(i,ncolumns-j-1) ] = data[ offset(i,j) ];
        }
    }
}

/** Reflect the contents of the left half to the right
    half. If n == nColumns(), then this sets

    matrix[i,j] = matrix[i,n-j]

    1 2 3      3 2 3
    4 5 6  =>  6 5 6
    7 8 9      9 8 9

    \throw SireError::incompatible_error
*/
void NMatrix::reflectRightToLeft()
{
    double *data = array.data();

    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncolumns/2; ++j)
        {
            data[ offset(i,j) ] = data[ offset(i,ncolumns-j-1) ];
        }
    }
}

/** Copy the contents of the top left diagonal to the bottom
    right diagonal. If n == nRows(), then this sets

    matrix[n-j,n-i] = matrix[i,j]

    1 2 3      1 2 3       [2,1] => [1,0]    [0,2] == [0,2]
    4 5 6  =>  4 5 2       [1,2] => [0,1]    [1,1] == [1,1]
    7 8 9      7 4 1       [2,2] => [0,0]    [2,0] == [2,0]

    This must be a square matrix.

    \throw SireError::incompatible_error
*/
void NMatrix::reflectTopLeftToBottomRight()
{
    assertSquare();

    double *data = array.data();

    for (int i=0; i<nrows-1; ++i)
    {
        for (int j=0; j<nrows-i-1; ++j)
        {
            data[ offset(nrows-j-1,nrows-i-1) ] = data[ offset(i,j) ];
        }
    }
}

/** Copy the contents of the top right diagonal to the bottom
    left diagonal. This sets matrix[j,i] = matrix[i,j]

    This must be a square matrix.

    1 2 3      1 2 3
    4 5 6  =>  2 5 6
    7 8 9      3 6 9

    \throw SireError::incompatible_error
*/
void NMatrix::reflectTopRightToBottomLeft()
{
    assertSquare();

    double *data = array.data();

    for (int i=0; i<nrows-1; ++i)
    {
        for (int j=i+1; j<nrows; ++j)
        {
            data[ offset(j,i) ] = data[ offset(i,j) ];
        }
    }
}

/** Copy the contents of the bottom right diagonal to the top
    left diagonal. If n == nRows(), then this sets

    matrix[i,j] = matrix[n-j,n-i]

    This must be a square matrix.

    1 2 3      9 6 3
    4 5 6  =>  8 5 6
    7 8 9      7 8 9

    \throw SireError::incompatible_error
*/
void NMatrix::reflectBottomRightToTopLeft()
{
    assertSquare();

    double *data = array.data();

    for (int i=0; i<nrows-1; ++i)
    {
        for (int j=0; j<nrows-i-1; ++j)
        {
            data[ offset(i,j) ] = data[ offset(nrows-j-1,nrows-i-1) ];
        }
    }
}

/** Copy the contents of the bottom left diagonal to the top
    right diagonal. This sets matrix[i,j] = matrix[j,i]

    This must be a square matrix.

    1 2 3      1 4 7
    4 5 6  =>  4 5 8
    7 8 9      7 8 9

    \throw SireError::incompatible_error
*/
void NMatrix::reflectBottomLeftToTopRight()
{
    assertSquare();

    double *data = array.data();

    for (int i=0; i<nrows-1; ++i)
    {
        for (int j=i+1; j<nrows; ++j)
        {
            data[ offset(i,j) ] = data[ offset(j,i) ];
        }
    }
}

/** Return the determinant of this matrix

    This uses LAPACK under the hood, for speed

    \throw SireError::incompatible_error
*/
double NMatrix::determinant() const
{
    this->assertSquare();

    std::pair< NMatrix,QVector<int> > factors_with_pivot = dgeco(*this);
    return dgedi_determinant(factors_with_pivot.first,
                             factors_with_pivot.second);
}

/** Return the trace of this matrix - this is only valid for a square matrix

    \throw SireError::incompatible_error
*/
double NMatrix::trace() const
{
    this->assertSquare();

    const double *d = array.constData();
    double sum = 0;

    for (int i=0; i<nrows; ++i)
    {
        sum += d[i*nrows + i];
    }

    return sum;
}

/** Return a vector containing the diagonal of this matrix - this is only
    valid for a square matrix

    \throw SireError::incompatible_error
*/
NVector NMatrix::diagonal() const
{
    this->assertSquare();

    const double *d = array.constData();

    NVector vector(nrows);
    double *v = vector.data();

    for (int i=0; i<nrows; ++i)
    {
        v[i] = d[i*nrows + i];
    }

    return vector;
}

/** Return whether or not this is a transposed matrix (data
    is stored in row-major order rather than column-major order) */
bool NMatrix::isTransposed() const
{
    return is_transpose;
}

/** Return the eigenvalues and eigenvectors of this matrix. This
    uses LAPACK under the hood for speed

    \throw SireError::incompatible_error
    \throw SireMaths::domain_error
*/
std::pair<NVector,NMatrix> NMatrix::diagonalise() const
{
    this->assertSquare();
    return dsyev(*this);
}
