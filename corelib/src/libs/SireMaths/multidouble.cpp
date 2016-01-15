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

#include "multidouble.h"

#include <QStringList>

#include "SireError/errors.h"

#include <QDebug>

using namespace SireMaths;

#ifdef MULTIFLOAT_AVX_IS_AVAILABLE
    static inline bool isAligned32(const void *pointer)
    {
        return (quintptr)pointer % size_t(32) == 0;
    }

    static void assertAligned32(const void *pointer, QString place)
    {
        if (not isAligned32(pointer))
            throw SireError::program_bug( QObject::tr(
                    "An unaligned MultiDouble has been created! %1")
                        .arg((quintptr)pointer % size_t(32)), place );
    }
#else
#ifdef MULTIFLOAT_SSE_IS_AVAILABLE
    static inline bool isAligned16(const void *pointer)
    {
        return (quintptr)pointer % size_t(16) == 0;
    }

    static void assertAligned16(const void *pointer, QString place)
    {
        if (not isAligned16(pointer))
            throw SireError::program_bug( QObject::tr(
                    "An unaligned MultiDouble has been created! %1")
                        .arg((quintptr)pointer % size_t(16)), place );
    }
#else
    static inline bool isAligned32(const void *pointer)
    {
        return (quintptr)pointer % size_t(32) == 0;
    }

    static void assertAligned32(const void *pointer, QString place)
    {
        if (not isAligned32(pointer))
            throw SireError::program_bug( QObject::tr(
                    "An unaligned MultiDouble has been created! %1")
                        .arg((quintptr)pointer % size_t(32)), place );
    }
#endif
#endif

void MultiDouble::assertAligned(const void *ptr, size_t size)
{
    if ( (quintptr)ptr % size != 0 )
        throw SireError::program_bug( QObject::tr(
                "An unaligned MultiDouble has been created! %1, %2, %3")
                    .arg((quintptr)ptr)
                    .arg((quintptr)ptr % size)
                    .arg(size), CODELOC );
}

/** Construct from the passed array. If size is greater than MultiDouble::size()
    then an error will be raised. If size is less than MultiDouble::size() then
    this vector will be padded with zeroes */
MultiDouble::MultiDouble(const double *array, int size)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        assertAligned32(this, CODELOC);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        assertAligned16(this, CODELOC);
    #else
        assertAligned32(this, CODELOC);
    #endif
    #endif

    if (size > MULTIFLOAT_SIZE)
        throw SireError::unsupported( QObject::tr(
                "Cannot fit an array of size %1 in this MultiDouble, as it is only "
                "capable of holding %2 values...").arg(size).arg(MULTIFLOAT_SIZE), CODELOC );

    if (size <= 0)
    {
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            v.x[0] = _mm256_set1_pd(0);
            v.x[1] = _mm256_set1_pd(0);
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            v.x[0] = _mm_set1_pd(0);
            v.x[1] = _mm_set1_pd(0);
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = 0;
            }
        #endif
        #endif
    }
    else if (size == MULTIFLOAT_SIZE)
    {
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            v.x[0] = _mm256_set_pd(array[3], array[2], array[1], array[0]);
            v.x[1] = _mm256_set_pd(array[7], array[6], array[5], array[4]);
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            //note that SSE packs things the 'wrong' way around
            v.x[0] = _mm_set_pd(array[1], array[0]);
            v.x[1] = _mm_set_pd(array[3], array[2]);
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = array[i];
            }
        #endif
        #endif
    }
    else
    {
        double tmp[MULTIFLOAT_SIZE];
        
        for (int i=0; i<size; ++i)
        {
            tmp[i] = array[i];
        }
        
        for (int i=size; i<MULTIFLOAT_SIZE; ++i)
        {
            tmp[i] = 0;
        }
        
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            v.x[0] = _mm256_set_pd(tmp[3], tmp[2], tmp[1], tmp[0]);
            v.x[1] = _mm256_set_pd(tmp[7], tmp[6], tmp[5], tmp[4]);
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            //note that sse packs things the 'wrong' way around
            v.x[0] = _mm_set_pd(tmp[1], tmp[0]);
            v.x[1] = _mm_set_pd(tmp[3], tmp[2]);
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
MultiDouble::MultiDouble(const QVector<float> &array)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        assertAligned32(this, CODELOC);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        assertAligned16(this, CODELOC);
    #else
        assertAligned32(this, CODELOC);
    #endif
    #endif

    QVector<double> darray;
    darray.reserve(array.count());
    
    for (int i=0; i<array.count(); ++i)
    {
        darray.append(array.constData()[i]);
    }

    this->operator=( MultiDouble(darray) );
}

/** Construct from the passed array - this must be the same size as the vector */
MultiDouble::MultiDouble(const QVector<double> &array)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        assertAligned32(this, CODELOC);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        assertAligned16(this, CODELOC);
    #else
        assertAligned32(this, CODELOC);
    #endif
    #endif

    this->operator=( MultiDouble(array.constData(), array.size()) );
}

/** Return whether or not this MultiDouble is correctly aligned. If it is not,
    then any SSE/AVX operations will fail */
bool MultiDouble::isAligned() const
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

QVector<MultiDouble> MultiDouble::fromArray(const double *array, int size)
{
    if (size == 0)
        return QVector<MultiDouble>();
    
    int nvecs = size / MULTIFLOAT_SIZE;
    int nremain = size % MULTIFLOAT_SIZE;
    
    QVector<MultiDouble> marray(nvecs + ( (nremain > 0) ? 1 : 0 ));
    
    MultiDouble *ma = marray.data();
    
    int idx = 0;
    
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        if (isAligned16(array))
        {
            for (int i=0; i<nvecs; ++i)
            {
                ma[i] = MultiDouble(array+idx, MULTIFLOAT_SIZE);
                idx += MULTIFLOAT_SIZE;
            }
    
            if (nremain > 0)
            {
                ma[marray.count()-1] = MultiDouble(array+idx, nremain);
            }
        }
        else
        {
            double _ALIGNED(16) tmp[MULTIFLOAT_SIZE];

            for (int i=0; i<nvecs; ++i)
            {
                for (int j=0; j<MULTIFLOAT_SIZE; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
            
                ma[i] = MultiDouble((double*)(&tmp), MULTIFLOAT_SIZE);
            }
            
            if (nremain > 0)
            {
                for (int j=0; j<nremain; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
                
                ma[marray.count()-1] = MultiDouble((double*)(&tmp), nremain);
            }
        }
    #else
        if (isAligned32(array))
        {
            for (int i=0; i<nvecs; ++i)
            {
                ma[i] = MultiDouble(array+idx, MULTIFLOAT_SIZE);
                idx += MULTIFLOAT_SIZE;
            }
    
            if (nremain > 0)
            {
                ma[marray.count()-1] = MultiDouble(array+idx, nremain);
            }
        }
        else
        {
            double _ALIGNED(32) tmp[MULTIFLOAT_SIZE];

            for (int i=0; i<nvecs; ++i)
            {
                for (int j=0; j<MULTIFLOAT_SIZE; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
            
                ma[i] = MultiDouble((double*)(&tmp), MULTIFLOAT_SIZE);
            }
            
            if (nremain > 0)
            {
                for (int j=0; j<nremain; ++j)
                {
                    tmp[j] = array[idx];
                    ++idx;
                }
                
                ma[marray.count()-1] = MultiDouble((double*)(&tmp), nremain);
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

/** Create an array of MultiFloats from the passed array of doubles. This
    will pad the end of the array with zeroes if necessary */
QVector<MultiDouble> MultiDouble::fromArray(const QVector<double> &array)
{
    return MultiDouble::fromArray(array.constData(), array.count());
}

QVector<MultiDouble> MultiDouble::fromArray(const float *array, int size)
{
    if (size == 0)
        return QVector<MultiDouble>();

    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        double _ALIGNED(16) tmp[MULTIFLOAT_SIZE];
    #else
        double _ALIGNED(32) tmp[MULTIFLOAT_SIZE];
    #endif
    
    int nvecs = size / MULTIFLOAT_SIZE;
    int nremain = size % MULTIFLOAT_SIZE;

    QVector<MultiDouble> marray(nvecs + ( (nremain > 0) ? 1 : 0 ));
    MultiDouble *a = marray.data();
    
    int idx = 0;
    
    for (int i=0; i<nvecs; ++i)
    {
        for (int j=0; j<MULTIFLOAT_SIZE; ++j)
        {
            tmp[j] = array[idx];
            ++idx;
        }
    
        a[i] = MultiDouble((double*)(&tmp), MULTIFLOAT_SIZE);
    }
    
    if (nremain > 0)
    {
        for (int j=0; j<nremain; ++j)
        {
            tmp[j] = array[idx];
            ++idx;
        }
        
        a[marray.count()-1] = MultiDouble((double*)(&tmp), nremain);
    }
    
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

/** Create an array of MultiFloats from the passed array of floats. This will
    pad the end of the array with zeroes if necessary */
QVector<MultiDouble> MultiDouble::fromArray(const QVector<float> &array)
{
    return MultiDouble::fromArray(array.constData(), array.count());
}

/** Return the passed MultiDouble converted back into a normal array */
QVector<double> MultiDouble::toArray(const QVector<MultiDouble> &array)
{
    if (array.isEmpty())
        return QVector<double>();
    
    QVector<double> ret;
    ret.reserve( array.count() * MULTIFLOAT_SIZE );
    
    for (int i=0; i<array.count(); ++i)
    {
        const MultiDouble &f = array.constData()[i];
        
        for (int j=0; j<MULTIFLOAT_SIZE; ++j)
        {
            ret.append(f[j]);
        }
    }
    
    return ret;
}

/** Return the passed MultiFloat converted back into a normal array of doubles */
QVector<double> MultiDouble::toDoubleArray(const QVector<MultiDouble> &array)
{
    return MultiDouble::toArray(array);
}

/** Comparison operator - only returns true if all elements are equal */
bool MultiDouble::operator==(const MultiDouble &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] != other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are not equal */
bool MultiDouble::operator!=(const MultiDouble &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] == other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are less */
bool MultiDouble::operator<(const MultiDouble &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] >= other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are greater */
bool MultiDouble::operator>(const MultiDouble &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] <= other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are less or equal */
bool MultiDouble::operator<=(const MultiDouble &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] > other.v.a[i])
            return false;
    }

    return true;
}

/** Comparison operator - only returns true if all elements are greater or equal */
bool MultiDouble::operator>=(const MultiDouble &other) const
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        if (v.a[i] < other.v.a[i])
            return false;
    }

    return true;
}

/** Return the ith value in the multifloat */
double MultiDouble::at(int i) const
{
    if (i < 0)
        i = MULTIFLOAT_SIZE + i;
    
    if (i < 0 or i >= MULTIFLOAT_SIZE)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot access element %1 of MultiDouble (holds only %2 values)")
                    .arg(i).arg(MULTIFLOAT_SIZE), CODELOC );
    }
    
    return v.a[i];
}

double MultiDouble::getitem(int i) const
{
    return at(i);
}

/** Negative operator */
MultiDouble MultiDouble::operator-() const
{
    MultiDouble ret;
    
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = -v.a[i];
    }
    
    return ret;
}

/** Set the ith value of the multifloat to 'value' */
void MultiDouble::set(int i, double value)
{
    if (i < 0)
        i = MULTIFLOAT_SIZE + i;
    
    if (i < 0 or i >= MULTIFLOAT_SIZE)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot access element %1 of MultiDouble (holds only %2 values)")
                    .arg(i).arg(MULTIFLOAT_SIZE), CODELOC );
    }

    v.a[i] = value;
}

/** Set the ith value without checking that i is valid */
void MultiDouble::quickSet(int i, double value)
{
    v.a[i] = value;
}

/** Return the ith value in the multifloat */
double MultiDouble::get(int i) const
{
    return at(i);
}

const char* MultiDouble::what() const
{
    return MultiDouble::typeName();
}

const char* MultiDouble::typeName()
{
    return "SireMaths::MultiDouble";
}

QString MultiDouble::toString() const
{
    QStringList vals;
    
    for (int i=0; i<this->count(); ++i)
    {
        vals.append( QString::number(v.a[i]) );
    }
    
    return QObject::tr("{ %1 }").arg(vals.join(", "));
}

QString MultiDouble::toBinaryString() const
{
    QStringList vals;
    
    for (int i=0; i<this->count(); ++i)
    {
        const unsigned char *c = reinterpret_cast<const unsigned char*>(&(v.a[i]));
        
        QString val("0x");
        
        for (unsigned int j=0; j<sizeof(double); ++j)
        {
            val.append( QString("%1").arg((unsigned short)(c[j]), 2, 16, QChar('0')) );
        }
        
        vals.append(val);
    }
    
    return QObject::tr("{ %1 }").arg(vals.join(", "));
}
