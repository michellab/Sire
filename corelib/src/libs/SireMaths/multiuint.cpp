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

#include "multiuint.h"

#include <QStringList>

#include "SireError/errors.h"

#include <QDebug>

using namespace SireMaths;

#ifdef MULTIFLOAT_AVX_IS_AVAILABLE
    static SIRE_ALWAYS_INLINE bool isAligned32(const void *pointer)
    {
        return (quintptr)pointer % size_t(32) == 0;
    }

    static void assertAligned32(const void *pointer, QString place)
    {
        if (not isAligned32(pointer))
            throw SireError::program_bug( QObject::tr(
                    "An unaligned MultiUInt has been created! %1")
                        .arg((quintptr)pointer % size_t(32)), place );
    }
#else
#ifdef MULTIFLOAT_SSE_IS_AVAILABLE
    static SIRE_ALWAYS_INLINE bool isAligned16(const void *pointer)
    {
        return (quintptr)pointer % size_t(16) == 0;
    }

    static void assertAligned16(const void *pointer, QString place)
    {
        if (not isAligned16(pointer))
            throw SireError::program_bug( QObject::tr(
                    "An unaligned MultiUInt has been created! %1")
                        .arg((quintptr)pointer % size_t(16)), place );
    }
#else
    static SIRE_ALWAYS_INLINE bool isAligned32(const void *pointer)
    {
        return (quintptr)pointer % size_t(32) == 0;
    }

    static void assertAligned32(const void *pointer, QString place)
    {
        if (not isAligned32(pointer))
            throw SireError::program_bug( QObject::tr(
                    "An unaligned MultiUInt has been created! %1")
                        .arg((quintptr)pointer % size_t(32)), place );
    }
#endif
#endif

void MultiUInt::assertAligned(const void *ptr, size_t size)
{
    if ( (quintptr)ptr % size != 0 )
        throw SireError::program_bug( QObject::tr(
                "An unaligned MultiUInt has been created! %1, %2, %3")
                    .arg((quintptr)ptr)
                    .arg((quintptr)ptr % size)
                    .arg(size), CODELOC );
}

/** Construct from the passed array. If size is greater than MultiInt::size()
    then an error will be raised. If size is less than MultiInt::size() then
    this float will be padded with zeroes */
MultiUInt::MultiUInt(const quint32 *array, int size)
{
    assertAligned();

    if (size > MULTIFLOAT_SIZE)
        throw SireError::unsupported( QObject::tr(
                "Cannot fit an array of size %1 in this MultiInt, as it is only "
                "capable of holding %2 values...").arg(size).arg(MULTIFLOAT_SIZE), CODELOC );

    if (size <= 0)
    {
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
                v.x = _mm256_set1_epi32(0);
            #else
                v.x[0] = _mm_set1_epi32(0);
                v.x[1] = _mm_set1_epi32(0);
            #endif
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            v.x = _mm_set1_epi32(0);
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = 0;
            }
        #endif
        #endif
    }
    else
    {
        quint32 tmp[MULTIFLOAT_SIZE];
        
        for (int i=0; i<size; ++i)
        {
            tmp[i] = *(reinterpret_cast<const qint32*>(&array[i]));
        }
        
        for (int i=size; i<MULTIFLOAT_SIZE; ++i)
        {
            tmp[i] = 0;
        }
        
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
                v.x = _mm256_set_epi32(tmp[7], tmp[6], tmp[5], tmp[4],
                                       tmp[3], tmp[2], tmp[1], tmp[0]);
            #else
                v.x[1] = _mm_set_epi32(tmp[7], tmp[6], tmp[5], tmp[4]);
                v.x[0] = _mm_set_epi32(tmp[3], tmp[2], tmp[1], tmp[0]);
            #endif
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            //note that sse packs things the 'wrong' way around
            v.x = _mm_set_epi32(tmp[3], tmp[2], tmp[1], tmp[0]);
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = tmp[i];
            }
        #endif
        #endif
    }
}

/** Construct from the passed array - this must be the same size as the vector */
MultiUInt::MultiUInt(const QVector<quint32> &array)
{
    assertAligned();
    this->operator=( MultiUInt(array.constData(), array.size()) );
}

/** Return whether or not this MultiInt is correctly aligned. If it is not,
    then any SSE operations will fail */
bool MultiUInt::isAligned() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return isAligned32(this);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return isAligned16(this);
    #else
        return true;
    #endif
    #endif
}

QVector<MultiUInt> MultiUInt::fromArray(const quint32 *array, int size)
{
    if (size == 0)
        return QVector<MultiUInt>();
    
    int nvecs = size / MULTIFLOAT_SIZE;
    int nremain = size % MULTIFLOAT_SIZE;
    
    QVector<MultiUInt> marray(nvecs + ( (nremain > 0) ? 1 : 0 ));
    MultiUInt *ma = marray.data();
    
    int idx = 0;
    
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        if (isAligned16(array))
        {
            for (int i=0; i<nvecs; ++i)
            {
                ma[i] = MultiUInt(array+idx, MULTIFLOAT_SIZE);
                idx += MULTIFLOAT_SIZE;
            }
    
            if (nremain > 0)
            {
                ma[marray.count()-1] = MultiUInt(array+idx, nremain);
            }
        }
        else
        {
            quint32 _ALIGNED(16) tmp[MULTIFLOAT_SIZE];

            for (int i=0; i<nvecs; ++i)
            {
                for (int j=0; j<MULTIFLOAT_SIZE; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
            
                ma[i] = MultiUInt((quint32*)(&tmp), MULTIFLOAT_SIZE);
            }
            
            if (nremain > 0)
            {
                for (int j=0; j<nremain; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
                
                ma[marray.count()-1] = MultiUInt((quint32*)(&tmp), nremain);
            }
        }
    #else
        if (isAligned32(array))
        {
            for (int i=0; i<nvecs; ++i)
            {
                ma[i] = MultiUInt(array+idx, MULTIFLOAT_SIZE);
                idx += MULTIFLOAT_SIZE;
            }
    
            if (nremain > 0)
            {
                ma[marray.count()-1] = MultiUInt(array+idx, nremain);
            }
        }
        else
        {
            quint32 _ALIGNED(32) tmp[MULTIFLOAT_SIZE];

            for (int i=0; i<nvecs; ++i)
            {
                for (int j=0; j<MULTIFLOAT_SIZE; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
            
                ma[i] = MultiUInt((quint32*)(&tmp), MULTIFLOAT_SIZE);
            }
            
            if (nremain > 0)
            {
                for (int j=0; j<nremain; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
                
                ma[marray.count()-1] = MultiUInt((quint32*)(&tmp), nremain);
            }
        }
    #endif

    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        assertAligned32(marray.constData(), CODELOC);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        assertAligned16(marray.constData(), CODELOC);
    #else
        assertAligned32(marray.constData(), CODELOC);
    #endif
    #endif
    
    return marray;
}

/** Create an array of MultiInts from the passed array of integers. This will
    pad the end of the array with zeroes if necessary */
QVector<MultiUInt> MultiUInt::fromArray(const QVector<quint32> &array)
{
    return MultiUInt::fromArray(array.constData(), array.count());
}

/** Return the passed MultiUInt converted back into a normal array */
QVector<quint32> MultiUInt::toArray(const QVector<MultiUInt> &array)
{
    if (array.isEmpty())
        return QVector<quint32>();
    
    QVector<quint32> ret;
    ret.reserve( array.count() * MULTIFLOAT_SIZE );
    
    for (int i=0; i<array.count(); ++i)
    {
        const MultiUInt &f = array.constData()[i];
        
        for (int j=0; j<MULTIFLOAT_SIZE; ++j)
        {
            ret.append(f[j]);
        }
    }
    
    return ret;
}

/** Comparison operator - only returns true if all elements are equal */
bool MultiUInt::operator==(const MultiUInt &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] != other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are not equal */
bool MultiUInt::operator!=(const MultiUInt &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] == other.v.a[i])
            return false;
    }

    return true;
}

/** Return whether all of the elements of this MultiInt are
    equal to 0x00000000 (e.g. every bit in the entire vector is 0) */
bool MultiUInt::isBinaryZero() const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        static const quint32 bin_zero = 0x00000000;
    
        if (*(reinterpret_cast<const quint32*>(&(v.a[i]))) != bin_zero)
            return false;
    }
    
    return true;
}

/** Return whether all of the elements of this MultiInt are
    not equal to 0x00000000 (e.g. at least one bit in the entire vector is 1) */
bool MultiUInt::isNotBinaryZero() const
{
    return not isBinaryZero();
}

/** Return whether or not at least one of the elements of this vector
    is binary zero (the float is equal to 0x00000000) */
bool MultiUInt::hasBinaryZero() const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        static const quint32 bin_zero = 0x00000000;
    
        if (*(reinterpret_cast<const quint32*>(&(v.a[i]))) == bin_zero)
            return true;
    }
    
    return false;
}

/** Return whether all of the elements of this MultiInt are
    equal to 0xFFFFFFFF (e.g. every bit in the entire vector is 1) */
bool MultiUInt::isBinaryOne() const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        static const quint32 bin_one = 0xFFFFFFFF;
    
        if (*(reinterpret_cast<const quint32*>(&(v.a[i]))) != bin_one)
            return false;
    }
    
    return true;
}

/** Return whether all of the elements of this MultiInt are
    not equal to 0xFFFFFFFF (e.g. at least one bit in the entire vector is 0) */
bool MultiUInt::isNotBinaryOne() const
{
    return not isBinaryOne();
}

/** Return whether or not at least one of the elements of this vector
    is binary one (the float is equal to 0xFFFFFFFF) */
bool MultiUInt::hasBinaryOne() const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        static const quint32 bin_one = 0xFFFFFFFF;
    
        if (*(reinterpret_cast<const quint32*>(&(v.a[i]))) == bin_one)
            return true;
    }
    
    return false;
}

/** Comparison operator - only returns true if all elements are less */
bool MultiUInt::operator<(const MultiUInt &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] >= other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are greater */
bool MultiUInt::operator>(const MultiUInt &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] <= other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are less or equal */
bool MultiUInt::operator<=(const MultiUInt &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] > other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are greater or equal */
bool MultiUInt::operator>=(const MultiUInt &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] < other.v.a[i])
            return false;
    }

    return true;
}

/** Return the ith value in the multifloat */
quint32 MultiUInt::at(int i) const
{
    if (i < 0)
        i = MULTIFLOAT_SIZE + i;
    
    if (i < 0 or i >= MULTIFLOAT_SIZE)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot access element %1 of MultiInt (holds only %2 values)")
                    .arg(i).arg(MULTIFLOAT_SIZE), CODELOC );
    }
    
    return v.a[i];
}

quint32 MultiUInt::getitem(int i) const
{
    return at(i);
}

/** Set the ith value of the MultiInt to 'value' */
void MultiUInt::set(int i, quint32 value)
{
    if (i < 0)
        i = MULTIFLOAT_SIZE + i;
    
    if (i < 0 or i >= MULTIFLOAT_SIZE)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot access element %1 of MultiInt (holds only %2 values)")
                    .arg(i).arg(MULTIFLOAT_SIZE), CODELOC );
    }

    v.a[i] = value;
}

/** Return the 
ith value in the MultiInt */
quint32 MultiUInt::get(int i) const
{
    return at(i);
}

const char* MultiUInt::what() const
{
    return MultiUInt::typeName();
}

const char* MultiUInt::typeName()
{
    return "SireMaths::MultiUInt";
}

QString MultiUInt::toString() const
{
    QStringList vals;
    
    for (int i=0; i<this->count(); ++i)
    {
        vals.append( QString::number(v.a[i]) );
    }
    
    return QObject::tr("{ %1 }").arg(vals.join(", "));
}

QString MultiUInt::toBinaryString() const
{
    QStringList vals;
    
    for (int i=0; i<this->count(); ++i)
    {
        const unsigned char *c = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        
        QString val("0x");
        
        for (unsigned int j=0; j<sizeof(qint32); ++j)
        {
            val.append( QString("%1").arg((unsigned short)(c[j]), 2, 16, QChar('0')) );
        }
        
        vals.append(val);
    }
    
    return QObject::tr("{ %1 }").arg(vals.join(", "));
}
