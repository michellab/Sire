/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "SireMaths/multifixed.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireMaths;

#define DOUBLE64 double

const DOUBLE64 double_magic = double(6755399441055744.0);  // 2^(52-16) * 1.5
                                                           // as double has 52 bits of mantissa

static QString toBinary(qint64 value)
{
    QStringList vals;
    
    const unsigned char *c = reinterpret_cast<const unsigned char*>(&(value));
        
    QString val("0x");
        
    for (unsigned int j=0; j<sizeof(qint64); ++j)
    {
        val.append( QString("%1").arg((unsigned short)(c[j]), 2, 16, QChar('0')) );
    }
    
    return val;
}

static const double SCALE_TO_FIXED = 2LL << 32;
static const double SCALE_FROM_FIXED = double(1) / SCALE_TO_FIXED;

static qint64 convertToFixed(DOUBLE64 value)
{
    // add the magic number to 'value'. This lines up the
    // decimal points so that the mantissa (lower 52bits of the double)
    // should contain the integer
    return qint64( (value * SCALE_TO_FIXED) + 0.5 );
}

static DOUBLE64 convertFromFixed(qint64 value)
{
    return double(value) * SCALE_FROM_FIXED;
}

/** Constructor - all of the elements are set to zero */
MultiFixed::MultiFixed()
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] = 0;
    }
}

/** Construct such that all elements are equal to 'value' */
MultiFixed::MultiFixed(double value)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] = convertToFixed(value);
    }
}

/** Construct from the passed array. If size is greater than MultiFixed::size()
    then an error will be raised. If size is less than MultiFixed::size() then
    this vector will be padded with zeroes */
MultiFixed::MultiFixed(const double *array, int size)
{
    if (size > MULTIFLOAT_SIZE)
        throw SireError::unsupported( QObject::tr(
                "Cannot fit an array of size %1 in this MultiFixed, as it is only "
                "capable of holding %2 values...").arg(size).arg(MULTIFLOAT_SIZE), CODELOC );

    if (size <= 0)
    {
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = 0;
        }
    }
    else if (size == MULTIFLOAT_SIZE)
    {
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = convertToFixed(array[i]);
        }
    }
    else
    {
        for (int i=0; i<size; ++i)
        {
            v.a[i] = convertToFixed(array[i]);
        }
        
        for (int i=size; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = 0;
        }
    }
}

/** Construct from the passed array. If size is greater than MultiFixed::size()
    then an error will be raised. If size is less than MultiFixed::size() then
    this vector will be padded with zeroes */
MultiFixed::MultiFixed(const QVector<double> &array)
{
    this->operator=( MultiFixed(array.constData(),array.count()) );
}

/** Construct from the passed MultiFloat */
MultiFixed::MultiFixed(const MultiFloat &value)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] = convertToFixed(value[i]);
    }
}

/** Copy constructor */
MultiFixed::MultiFixed(const MultiFixed &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] = other.v.a[i];
    }
}

/** Destructor */
MultiFixed::~MultiFixed()
{}

/** Convert the passed array of doubles to an array of MultiFixed values.
    Note that the array may be returned padded with zeroes */
QVector<MultiFixed> MultiFixed::fromArray(const QVector<double> &array)
{
    if (array.isEmpty())
        return QVector<MultiFixed>();

    QVector<MultiFixed> marray;
    
    int nvecs = array.count() / MULTIFLOAT_SIZE;
    int nremain = array.count() % MULTIFLOAT_SIZE;
    
    marray.reserve(nvecs + ( (nremain == 1) ? 1 : 0 ));
    
    int idx = 0;
    
    for (int i=0; i<nvecs; ++i)
    {
        marray.append( MultiFixed((double*)(&(array.constData()[idx])), MULTIFLOAT_SIZE) );
        idx += MULTIFLOAT_SIZE;
    }
    
    if (nremain > 0)
    {
        marray.append( MultiFixed((double*)(&(array.constData()[idx])), nremain) );
    }
    
    return marray;
}

/** Return convert the passed MultiFixed array back into an array of doubles */
QVector<double> MultiFixed::toArray(const QVector<MultiFixed> &array)
{
    if (array.isEmpty())
        return QVector<double>();
    
    QVector<double> ret;
    ret.reserve( array.count() * MULTIFLOAT_SIZE );
    
    for (int i=0; i<array.count(); ++i)
    {
        const MultiFixed &f = array.constData()[i];
        
        for (int j=0; j<MULTIFLOAT_SIZE; ++j)
        {
            ret.append(f[j]);
        }
    }
    
    return ret;
}

/** Copy assignment operator */
MultiFixed& MultiFixed::operator=(const MultiFixed &other)
{
    if (this != &other)
    {
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    }
    
    return *this;
}

/** Comparison operator. Returns true if all values in this vector are
    equal to the corresponding values in other */
bool MultiFixed::operator==(const MultiFixed &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] != other.v.a[i])
            return false;
    }
    
    return true;
}

/** Comparison operator */
bool MultiFixed::operator!=(const MultiFixed &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] == other.v.a[i])
            return false;
    }
    
    return true;
}

/** Less than operator */
bool MultiFixed::operator<(const MultiFixed &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] >= other.v.a[i])
            return false;
    }
    
    return true;
}

/** Greater than operator */
bool MultiFixed::operator>(const MultiFixed &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] <= other.v.a[i])
            return false;
    }
    
    return true;
}

/** Less or equal operator */
bool MultiFixed::operator<=(const MultiFixed &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] > other.v.a[i])
            return false;
    }
    
    return true;
}

/** Greater or equal operator */
bool MultiFixed::operator>=(const MultiFixed &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] < other.v.a[i])
            return false;
    }
    
    return true;
}

static qint64 getBinaryOne()
{
    const quint64 x = 0xFFFFFFFFFFFFFFFFULL;
    return *(reinterpret_cast<const qint64*>(&x));
}

#define MULTIFIXED_BINONE getBinaryOne()

/** Compare each element of the two vectors. Return 0x00000000000000000 if
    the element is not equal, 0x1111111111111111 if they are */
MultiFixed MultiFixed::compareEqual(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = (v.a[i] == other.v.a[i]) ? MULTIFIXED_BINONE : 0x0;
    }
    
    return ret;
}

/** Compare each element for inequality */
MultiFixed MultiFixed::compareNotEqual(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIFIXED_BINONE : 0x0;
    }
    
    return ret;
}

/** Compare each element for less */
MultiFixed MultiFixed::compareLess(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIFIXED_BINONE : 0x0;
    }
    
    return ret;
}

/** Compare each element for greater */
MultiFixed MultiFixed::compareGreater(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIFIXED_BINONE : 0x0;
    }
    
    return ret;
}

/** Compare each element for less or equal */
MultiFixed MultiFixed::compareLessEqual(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIFIXED_BINONE : 0x0;
    }
    
    return ret;
}

/** Compare each element for greater or equal */
MultiFixed MultiFixed::compareGreaterEqual(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIFIXED_BINONE : 0x0;
    }
    
    return ret;
}

const char* MultiFixed::what() const
{
    return MultiFixed::typeName();
}

const char* MultiFixed::typeName()
{
    return "SireMaths::MultiFixed";
}

QString MultiFixed::toString() const
{
    QStringList vals;
    
    for (int i=0; i<this->count(); ++i)
    {
        vals.append( QString::number(this->get(i)) );
    }
    
    return QObject::tr("{ %1 }").arg(vals.join(", "));
}

QString MultiFixed::toBinaryString() const
{
    QStringList vals;
    
    for (int i=0; i<this->count(); ++i)
    {
        const unsigned char *c = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        
        QString val("0x");
        
        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            val.append( QString("%1").arg((unsigned short)(c[j]), 2, 16, QChar('0')) );
        }
        
        vals.append(val);
    }
    
    return QObject::tr("{ %1 }").arg(vals.join(", "));
}

/** Return the number of elements in the vector */
int MultiFixed::size()
{
    return MULTIFLOAT_SIZE;
}

/** Return the number of elements in the vector */
int MultiFixed::count()
{
    return MULTIFLOAT_SIZE;
}

/** Return the value of the ith element of the vector */
double MultiFixed::operator[](int i) const
{
    if (i < 0)
        i = MULTIFLOAT_SIZE + i;
    
    if (i < 0 or i >= MULTIFLOAT_SIZE)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot access element %1 of MultiFixed (holds only %2 values)")
                    .arg(i).arg(MULTIFLOAT_SIZE), CODELOC );
    }
    
    return convertFromFixed(v.a[i]);
}

/** Set the ith element of this vector to 'value' */
void MultiFixed::set(int i, double value)
{
    if (i < 0)
        i = MULTIFLOAT_SIZE + i;
    
    if (i < 0 or i >= MULTIFLOAT_SIZE)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot access element %1 of MultiFixed (holds only %2 values)")
                    .arg(i).arg(MULTIFLOAT_SIZE), CODELOC );
    }

    v.a[i] = convertToFixed(value);
}

/** Return the value of the ith element of the vector */
double MultiFixed::get(int i) const
{
    return this->operator[](i);
}

/** Negate this vector */
MultiFixed MultiFixed::operator-() const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = -v.a[i];
    }
    
    return ret;
}

/** Addition operator */
MultiFixed MultiFixed::operator+(const MultiFixed &other) const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = v.a[i] + other.v.a[i];
    }
    
    return ret;
}

/** Subtraction operator */
MultiFixed MultiFixed::operator-(const MultiFixed &other) const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = v.a[i] - other.v.a[i];
    }
    
    return ret;
}

/** Multiplication operator */
MultiFixed MultiFixed::operator*(const MultiFixed &other) const
{
    //due to the high chance of overflow it is probably easier
    //to do this by back converting to a double
    
    double ret[ MULTIFLOAT_SIZE ];
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret[i] = convertFromFixed(v.a[i]) * convertFromFixed(other.v.a[i]);
    }
    
    return MultiFixed( (double*)(&ret), MULTIFLOAT_SIZE );
}

/** Division operator */
MultiFixed MultiFixed::operator/(const MultiFixed &other) const
{
    //due to the high chance of overflow it is probably easier
    //to do this by back converting to a double
    
    double ret[ MULTIFLOAT_SIZE ];
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret[i] = convertFromFixed(v.a[i]) / convertFromFixed(other.v.a[i]);
    }
    
    return MultiFixed( (double*)(&ret), MULTIFLOAT_SIZE );
}

/** In-place addition operator */
MultiFixed& MultiFixed::operator+=(const MultiFixed &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] += other.v.a[i];
    }
    
    return *this;
}

/** In-place subtraction operator */
MultiFixed& MultiFixed::operator-=(const MultiFixed &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] -= other.v.a[i];
    }
    
    return *this;
}

/** In-place multiplication operator */
MultiFixed& MultiFixed::operator*=(const MultiFixed &other)
{
    this->operator=( this->operator*(other) );
    
    return *this;
}

/** In-place division operator */
MultiFixed& MultiFixed::operator/=(const MultiFixed &other)
{
    this->operator=( this->operator/(other) );
    return *this;
}

/** Logical bitwise not */
MultiFixed MultiFixed::operator!() const
{
    return this->logicalNot();
}

/** Logical bitwise and operator */
MultiFixed MultiFixed::operator&(const MultiFixed &other) const
{
    return this->logicalAnd(other);
}

/** Logical bitwise or operator */
MultiFixed MultiFixed::operator|(const MultiFixed &other) const
{
    return this->logicalOr(other);
}

/** Logical bitwise xor operator */
MultiFixed MultiFixed::operator^(const MultiFixed &other) const
{
    return this->logicalXor(other);
}

/** In-place logical bitwise and operator */
MultiFixed& MultiFixed::operator&=(const MultiFixed &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v
                    = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            char_v[j] &= other_char_v[j];
        }
    }
    
    return *this;
}

/** In-place logical bitwise or operator */
MultiFixed& MultiFixed::operator|=(const MultiFixed &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v
                    = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            char_v[j] |= other_char_v[j];
        }
    }
    
    return *this;
}

/** In-place logical bitwise xor operator */
MultiFixed& MultiFixed::operator^=(const MultiFixed &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v
                    = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            char_v[j] ^= other_char_v[j];
        }
    }
    
    return *this;
}

/** Logical bitwise not operator */
MultiFixed MultiFixed::logicalNot() const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            ret_char_v[j] = !char_v[j];
        }
    }

    return ret;
}

/** Logical bitwise and operator */
MultiFixed MultiFixed::logicalAnd(const MultiFixed &other) const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v =
                                    reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            ret_char_v[j] = char_v[j] & other_char_v[j];
        }
    }

    return ret;
}

/** Logical bitwise and not */
MultiFixed MultiFixed::logicalAndNot(const MultiFixed &other) const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v =
                                    reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            ret_char_v[j] = !(char_v[j] & other_char_v[j]);
        }
    }

    return ret;
}

/** Logical bitwise or */
MultiFixed MultiFixed::logicalOr(const MultiFixed &other) const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v =
                                    reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            ret_char_v[j] = char_v[j] | other_char_v[j];
        }
    }

    return ret;
}

/** Logical bitwise xor */
MultiFixed MultiFixed::logicalXor(const MultiFixed &other) const
{
    MultiFixed ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        const unsigned char *other_char_v =
                                    reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

        for (unsigned int j=0; j<sizeof(qint64); ++j)
        {
            ret_char_v[j] = char_v[j] ^ other_char_v[j];
        }
    }

    return ret;
}

/** Multiply 'val0' and 'val1' and add it onto this vector */
MultiFixed& MultiFixed::multiplyAdd(const MultiFixed &val0, const MultiFixed &val1)
{
    return this->operator+=( val0 * val1 );
}

/** Return the max of this vector with other */
MultiFixed MultiFixed::max(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = ( v.a[i] >= other.v.a[i] ? v.a[i] : other.v.a[i] );
    }
    
    return ret;
}

/** Return the min of this vector with other */
MultiFixed MultiFixed::min(const MultiFixed &other) const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = ( v.a[i] <= other.v.a[i] ? v.a[i] : other.v.a[i] );
    }
    
    return ret;
}

/** Return the reciprocal of this number */
MultiFixed MultiFixed::reciprocal() const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = convertToFixed( double(1) / convertFromFixed(v.a[i]) );
    }
    
    return ret;
}

/** Return the square root of this number */
MultiFixed MultiFixed::sqrt() const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = convertToFixed( std::sqrt(convertFromFixed(v.a[i])) );
    }
    
    return ret;
}

/** Return the reciprocal square root of this number */
MultiFixed MultiFixed::rsqrt() const
{
    MultiFixed ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = convertToFixed( double(1) / std::sqrt(convertFromFixed(v.a[i])) );
    }
    
    return ret;
}

/** Rotate this vector in the same direction as MultiFloat::rotate() */
MultiFixed MultiFixed::rotate() const
{
    MultiFixed ret;
    
    for (int i=1; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i-1] = v.a[i];
    }
    
    ret.v.a[MULTIFLOAT_SIZE-1] = v.a[0];
    
    return ret;
}

/** Return the sum of all of the elements of this vector */
double MultiFixed::sum() const
{
    //form this sum in blocks of four. This should help preserve
    //the order of addition in case the size of MultiFloat changes...
    double total = 0;
    
    for (int i=0; i<MULTIFLOAT_SIZE; i+=4)
    {
        total += ( convertFromFixed(v.a[i]) +
                   convertFromFixed(v.a[i+1]) +
                   convertFromFixed(v.a[i+2]) +
                   convertFromFixed(v.a[i+3]) );
    }
    
    return total;
}
