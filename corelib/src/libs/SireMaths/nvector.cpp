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

#include "nvector.h"
#include "nmatrix.h"
#include "vector.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "SireMaths/errors.h"

using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<NVector> r_nvector(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NVector &vector)
{
    writeHeader(ds, r_nvector, 1);
    
    SharedDataStream sds(ds);
    
    sds << vector.array;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NVector &vector)
{
    VersionID v = readHeader(ds, r_nvector);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> vector.array;
    }
    else
        throw version_error(v, "1", r_nvector, CODELOC);
        
    return ds;
}

/** Null constructor */
NVector::NVector()
{}

/** Construct a column vector with 'nrows' rows - the values
    of the rows are not initialised */
NVector::NVector(int nrows) : array(nrows)
{
    array.squeeze();
}

/** Construct a column vector with 'nrows' rows, all
    initialised to the value 'initial_value' */
NVector::NVector(int nrows, double initial_value)
        : array(nrows, initial_value)
{
    array.squeeze();
}

/** Construct from the passed 3-vector */
NVector::NVector(const Vector &vector) : array(3)
{
    array.squeeze();
    
    memcpy(array.data(), vector.constData(), 3*sizeof(double));
}

/** Construct from the passed vector */
NVector::NVector(const QVector<double> &vector) : array(vector)
{
    array.squeeze();
}

/** Copy constructor */
NVector::NVector(const NVector &other) : array(other.array)
{}

/** Destructor */
NVector::~NVector()
{}

const char* NVector::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NVector>() );
}

const char* NVector::what() const
{
    return NVector::typeName();
}

/** Copy assignment operator */
NVector& NVector::operator=(const NVector &other)
{
    array = other.array;
    return *this;
}

/** Comparison operator */
bool NVector::operator==(const NVector &other) const
{
    return array == other.array;
}

/** Comparison operator */
bool NVector::operator!=(const NVector &other) const
{
    return array != other.array;
}

/** Assert that the index 'i' is valid

    \throw SireError::invalid_index
*/
void NVector::assertValidIndex(int i) const
{
    if (i < 0 or i >= array.count())
        throw SireError::invalid_index( QObject::tr(
                "The column vector with dimension %1 does not have a value "
                "at index %2.")
                    .arg(array.count()).arg(i), CODELOC );
}

/** Assert that the index [i,j] is valid

    \throw SireError::invalid_index
*/
void NVector::assertValidIndex(int i, int j) const
{
    if (j != 0)
        throw SireError::invalid_index( QObject::tr(
                "A column vector only has a single column, so the vector "
                "with dimension [%1,1] cannot have an index [%2,%3].")
                    .arg(array.count()).arg(i).arg(j), CODELOC );

    if (i < 0 or i >= array.count())
        throw SireError::invalid_index( QObject::tr(
                "The column vector with dimension [%1,1] does not have a value "
                "at index [%2,%3].")
                    .arg(array.count()).arg(i).arg(j), CODELOC );
}

/** Assert that this column vector has 'nrows' rows

    \throw SireError::incompatible_error
*/
void NVector::assertNRows(int nrows) const
{
    if (nrows != array.count())
        throw SireError::incompatible_error( QObject::tr(
                "The column vector with dimension %1 is not compatible "
                "as the number of rows does not equal %2.")
                    .arg(array.count()).arg(nrows), CODELOC );
}

/** Assert that this column vector has 'ncolumns' columns - note
    that this is a column vector, so only has 1 column! 
    
    \throw SireError::incompatible_error
*/
void NVector::assertNColumns(int ncolumns) const
{
    if (array.isEmpty())
        throw SireError::incompatible_error( QObject::tr(
                "The null vector is incompatible as the number of columns "
                "does not equal %1.").arg(ncolumns), CODELOC );

    if (ncolumns != 1)
        throw SireError::incompatible_error( QObject::tr(
                "By definition, a column vector of dimension %1 only "
                "has one column, so it is incompatible with something that "
                "requires the number of columns to be equal to %2.")
                    .arg(array.count()).arg(ncolumns), CODELOC );
}

/** Return a modifiable reference to the ith row of this 
    column vector
    
    \throw SireError::invalid_index
*/
double& NVector::operator[](int i)
{
    this->assertValidIndex(i);
    return array.data()[i];
}

/** Return a modifiable reference to the ith row of this 
    column vector
    
    \throw SireError::invalid_index
*/
double& NVector::operator()(int i)
{
    return this->operator[](i);
}

/** Return a modifiable reference to the value at [i,j]
    
    \throw SireError::invalid_index
*/
double& NVector::operator()(int i, int j)
{
    this->assertValidIndex(i,j);
    return array.data()[i];
}

/** Return the value at index 'i'

    \throw SireError::invalid_index
*/
const double& NVector::operator[](int i) const
{
    this->assertValidIndex(i);
    return array.constData()[i];
}

/** Return the value at index 'i'

    \throw SireError::invalid_index
*/
const double& NVector::operator()(int i) const
{
    return this->operator[](i);
}

/** Return the value at index [i,j] 

    \throw SireError::invalid_index
*/
const double& NVector::operator()(int i, int j) const
{
    this->assertValidIndex(i,j);
    return array.constData()[i];
}

/** Return a raw pointer to the data of this vector */
double* NVector::data()
{
    return array.data();
}

/** Return a raw pointer to the data of this vector */
const double* NVector::data() const
{
    return array.constData();
}

/** Return a raw pointer to the data of this vector */
const double* NVector::constData() const
{
    return array.constData();
}

/** Return the number of rows in this column vector */
int NVector::nRows() const
{
    return array.count();
}

/** Return the number of columns in this column vector (0 or 1) */
int NVector::nColumns() const
{
    if (array.isEmpty())
        return 0;
    else
        return 1;
}

/** Vector addition */
NVector& NVector::operator+=(const NVector &other)
{
    this->assertNRows(other.nRows());
    
    double *d = array.data();
    const double *other_d = other.array.constData();
    
    const int sz = array.count();
    
    for (int i=0; i<sz; ++i)
    {
        d[i] += other_d[i];
    }
    
    return *this;
}

/** Vector subtraction */
NVector& NVector::operator-=(const NVector &other)
{
    this->assertNRows(other.nRows());
    
    double *d = array.data();
    const double *other_d = other.array.constData();
    
    const int sz = array.count();
    
    for (int i=0; i<sz; ++i)
    {
        d[i] -= other_d[i];
    }
    
    return *this;
}

/** Vector scale */
NVector& NVector::operator*=(double scale)
{
    double *d = array.data();
    const int sz = array.count();

    if (scale == 0)
    {
        for (int i=0; i<sz; ++i)
        {
            d[i] = 0;
        }
    }
    else
    {
        for (int i=0; i<sz; ++i)
        {
            d[i] *= scale;
        }
    }
    
    return *this;
}

/** Vector scale 

    \throw SireMaths::domain_error
*/
NVector& NVector::operator/=(double scale)
{
    if (scale == 0)
        throw SireMaths::domain_error( QObject::tr(
                "This code cannot handle dividing a vector by 0."), CODELOC );

    return this->operator*=( 1 / scale );
}

/** Return the negative of the vector */
NVector NVector::operator-() const
{
    NVector ret(*this);
    
    const int sz = array.count();
    double *d = ret.array.data();
    
    for (int i=0; i<sz; ++i)
    {
        d[i] = -d[i];
    }
    
    return ret;
}

/** Vector addition */
NVector NVector::operator+(const NVector &other) const
{
    NVector ret(*this);
    ret += other;
    
    return ret;
}

/** Vector subtraction */
NVector NVector::operator-(const NVector &other) const
{
    NVector ret(*this);
    ret -= other;
    return ret;
}

/** Vector scale */
NVector NVector::operator*(double scale) const
{
    NVector ret(*this);
    ret *= scale;
    return ret;
}

/** Vector scale */
NVector NVector::operator/(double scale) const
{
    NVector ret(*this);
    ret /= scale;
    return *this;
}

/** Set the value at [i,0] to 'value'

    \throw SireError::invalid_index
*/
void NVector::set(int i, double value)
{
    this->assertValidIndex(i);
    array.data()[i] = value;
}

/** Set the value at [i,j] to 'value'

    \throw SireError::invalid_index
*/
void NVector::set(int i, int j, double value)
{
    this->assertValidIndex(i,j);
    array.data()[i] = value;
}

/** Set all values in this vector to 'value' */
void NVector::setAll(double value)
{
    const int sz = array.count();
    double *d = array.data();
    
    for (int i=0; i<sz; ++i)
    {
        d[i] = value;
    }
}

/** Return the number of values in this vector */
int NVector::count() const
{
    return array.count();
}

/** Return the number of values in this vector */
int NVector::size() const
{
    return this->count();
}

/** Return the length squared of this vector */
double NVector::length2() const
{
    const int sz = array.count();
    const double *d = array.constData();
    
    double l2(0);
    
    for (int i=0; i<sz; ++i)
    {
        l2 += d[i]*d[i];
    }
    
    return l2;
}

/** Return the length of this vector */
double NVector::length() const
{
    return std::sqrt( this->length2() );
}

/** Normalise this vector

    \throw SireMaths::domain_error
*/
NVector NVector::normalise() const
{
    double l = this->length();
    
    if (l == 0)
        throw SireMaths::domain_error( QObject::tr(
                "It is not possible to normalise a vector of length 0."),
                    CODELOC );
                    
    l = 1 / l;
    
    return *this * l;
}

/** Return the sum of the elements of this vector */
double NVector::sum() const
{
    const double *d = array.constData();
    
    int sz = array.count();
    
    double total = 0;
    
    for (int i=0; i<sz; ++i)
    {
        total += d[i];
    }
    
    return total;
}

/** Return whether or not this is a zero vector */
bool NVector::isZero() const
{
    return this->length2() == 0;
}

/** Return a string representation of this vector */
QString NVector::toString() const
{
    if (array.count() == 0)
        return "( )";
        
    else if (array.count() == 1)
    {
        const double *d = array.constData();
    
        return QString("( %1 )").arg( QString::number(d[0], 'g', 8) );
    }

    QStringList rows;
    
    const double *d = array.constData();
    
    for (qint32 i=0; i<array.count(); ++i)
    {
        QString numstr = QString::number(d[i], 'g', 6);
    
        if (i == 0)
            rows.append( QString("/ %1 \\").arg(numstr, 8) );
        else if (i == this->nRows() - 1)
            rows.append( QString("\\ %1 /").arg(numstr, 8) );
        else
            rows.append( QString("| %1 |").arg(numstr, 8) );
    }
    
    return rows.join("\n");
}

/** Return the dot product of this vector with 'other'

    \throw SireError::incompatible_error
*/
double NVector::dot(const NVector &other) const
{
    this->assertNRows(other.nRows());
    
    const int sz = array.count();
    const double *d = array.constData();
    const double *other_d = other.array.constData();
    
    double sum(0);
    
    for (int i=0; i<sz; ++i)
    {
        sum += d[i]*other_d[i];
    }
    
    return sum;
}

/** Return the cross product of this vector with 'other'

    \throw SireError::incompatible_error
    \throw SireMaths::domain_error
*/
NVector NVector::cross(const NVector &other) const
{
    this->assertNRows(other.nRows());
    
    if (this->nRows() == 3)
    {
        const double *d = array.constData();
        const double *other_d = other.array.constData();
    
        return Vector::cross( Vector(d[0], d[1], d[2]),
                              Vector(other_d[0], other_d[1], other_d[2]) );
    }
    else if (this->nRows() == 7)
    {
        throw SireError::incomplete_code( QObject::tr(
                "While the cross product is defined for 7-dimensional vectors, "
                "it is not yet implemented in this code!"), CODELOC );
    }
    else
    {
        throw SireError::incompatible_error( QObject::tr(
                "The cross-product is only defined for 3- and 7-dimensional "
                "vectors. A column vector of rank %1 is not supported.")
                    .arg(this->nRows()), CODELOC );
    }
    
    return NVector();
}

/** Return the transpose of this column vector (a row vector, which
    is implemented in this code as a NMatrix) */
NMatrix NVector::transpose() const
{
    return NMatrix(*this, true);
}
