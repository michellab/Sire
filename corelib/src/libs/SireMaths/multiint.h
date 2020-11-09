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

#ifndef SIREMATHS_MULTIINT_H
#define SIREMATHS_MULTIINT_H

#include "multifloat.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

/** This class provides a vectorised 32bit signed integer. This represents
    a single vector of integers on the compiled machine, e.g.
    4 integers if we use SSE2, 8 integers for AVX/AVX2
    (note that AVX represents it as two SSE vectors, while AVX2 
     uses a single AVX vector)
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT MultiInt
{
public:
    MultiInt();
    
    MultiInt(qint32 value);
    
    MultiInt(const qint32 *array, int size);
    MultiInt(const QVector<qint32> &array);
    MultiInt(const MultiFloat &other);
    
    MultiInt(const MultiInt &other);
    
    ~MultiInt();
    
    bool isAligned() const;
    
    static QVector<MultiInt> fromArray(const QVector<qint32> &array);
    static QVector<MultiInt> fromArray(const qint32 *array, int size);
    
    static QVector<qint32> toArray(const QVector<MultiInt> &array);
    
    MultiInt& operator=(const MultiInt &other);
    MultiInt& operator=(qint32 value);
    MultiInt& operator=(const MultiFloat &other);
    
    bool operator==(const MultiInt &other) const;
    bool operator!=(const MultiInt &other) const;
    
    bool operator<(const MultiInt &other) const;
    bool operator>(const MultiInt &other) const;
    
    bool operator<=(const MultiInt &other) const;
    bool operator>=(const MultiInt &other) const;
    
    MultiInt compareEqual(const MultiInt &other) const;
    MultiInt compareNotEqual(const MultiInt &other) const;

    MultiInt compareLess(const MultiInt &other) const;
    MultiInt compareGreater(const MultiInt &other) const;
    
    MultiInt compareLessEqual(const MultiInt &other) const;
    MultiInt compareGreaterEqual(const MultiInt &other) const;
    
    MultiFloat reinterpretCastToFloat() const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    QString toBinaryString() const;
    
    static int size();
    static int count();
    
    qint32 operator[](int i) const;
    
    void set(int i, qint32 value);
    qint32 get(int i) const;
    
    qint32 at(int i) const;
    qint32 getitem(int i) const;
    
    MultiInt operator-() const;
    
    MultiInt operator+(const MultiInt &other) const;
    MultiInt operator-(const MultiInt &other) const;
    MultiInt operator*(const MultiInt &other) const;
    
    MultiInt& operator+=(const MultiInt &other);
    MultiInt& operator-=(const MultiInt &other);
    MultiInt& operator*=(const MultiInt &other);
    
    MultiInt operator!() const;
    MultiInt operator&(const MultiInt &other) const;
    MultiInt operator|(const MultiInt &other) const;
    MultiInt operator^(const MultiInt &other) const;

    MultiInt& operator&=(const MultiInt &other);
    MultiInt& operator|=(const MultiInt &other);
    MultiInt& operator^=(const MultiInt &other);

    MultiInt logicalNot() const;
    
    MultiInt logicalAnd(const MultiInt &other) const;
    MultiInt logicalAndNot(const MultiInt &other) const;
    
    MultiInt logicalOr(const MultiInt &other) const;
    MultiInt logicalXor(const MultiInt &other) const;
    
    MultiInt max(const MultiInt &other) const;
    MultiInt min(const MultiInt &other) const;
    
    MultiInt rotate() const;
    
    bool isBinaryZero() const;
    bool isNotBinaryZero() const;
    bool hasBinaryZero() const;
    
    bool isBinaryOne() const;
    bool isNotBinaryOne() const;
    bool hasBinaryOne() const;
    
    qint32 sum() const;
    qint64 doubleSum() const;

private:
    /* Make the other Multi?? classes friends, so that they
       can read the vector data directly */
    friend class MultiFloat;
    friend class MultiUInt;

    static void assertAligned(const void *ptr, size_t size);

    #ifndef SIRE_SKIP_INLINE_FUNCTIONS
        #ifndef MULTIFLOAT_CHECK_ALIGNMENT
            void assertAligned(){}
        #endif

        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
                _ALIGNED(32) union
                {
                    __m256i x;
                    qint32 a[8];
                } v;
    
                MultiInt(__m256i val)
                {
                    v.x = val;
                }
            #else
                _ALIGNED(32) union
                {
                    __m128i x[2];
                    qint32 a[8];
                } v;
        
                MultiInt(__m128i val0, __m128i val1)
                {
                    v.x[0] = val0;
                    v.x[1] = val1;
                }
            #endif

            #define MULTIINT_BINONE getBinaryOne()

            static qint32 getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return *(reinterpret_cast<const qint32*>(&x));
            }

            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if (quintptr(this) % 32 != 0)
                        assertAligned(this, 32);
                }
            #endif
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            _ALIGNED(16) union
            {
                __m128i x;
                qint32 a[4];
            } v;

            MultiInt(__m128i sse_val)
            {
                v.x = sse_val;
            }

            #define MULTIINT_BINONE getBinaryOne()

            static qint32 getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return *(reinterpret_cast<const qint32*>(&x));
            }

            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if (quintptr(this) % 16 != 0)
                        assertAligned(this, 16);
                }
            #endif
        #else
            _ALIGNED(32) union
            {
                qint32 a[MULTIFLOAT_SIZE];
            } v;
            #define MULTIINT_BINONE getBinaryOne()

            static qint32 getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return *(reinterpret_cast<const qint32*>(&x));
            }
            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if (quintptr(this) % 32 != 0)
                        assertAligned(this, 32);
                }
            #endif
        #endif
        #endif
    #endif

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor. This creates a MultiInt with an undefined initial state */
SIRE_ALWAYS_INLINE
MultiInt::MultiInt()
{
    assertAligned();

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

/** Construct a MultiFloat with all values equal to 'val' */
SIRE_ALWAYS_INLINE
MultiInt::MultiInt(qint32 val)
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_set1_epi32(val);
        #else
            v.x[0] = _mm_set1_epi32(val);
            v.x[1] = _mm_set1_epi32(val);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_set1_epi32(val);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = val;
        }
    #endif
    #endif
}

/** Copy constructor */
SIRE_ALWAYS_INLINE
MultiInt::MultiInt(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = other.v.x;
        #else
            v.x[0] = other.v.x[0];
            v.x[1] = other.v.x[1];
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = other.v.x;
    #else
       for (int i=0; i<MULTIFLOAT_SIZE; ++i)
       {
           v.a[i] = other.v.a[i];
       }
    #endif
    #endif
}

/** Return the ith value in the MultiInt - note that this is
    a quick function that does no bounds checking */
SIRE_ALWAYS_INLINE qint32 MultiInt::operator[](int i) const
{
    return v.a[i];
}

/** Assignment operator */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator=(const MultiInt &other)
{
    if (this != &other)
    {
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
                v.x = other.v.x;
            #else
                v.x[0] = other.v.x[0];
                v.x[1] = other.v.x[1];
            #endif
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            v.x = other.v.x;
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = other.v.a[i];
            }
        #endif
        #endif
    }
    
    return *this;
}

/** Assignment operator */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator=(qint32 value)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_set1_epi32(value);
        #else
            v.x[0] = _mm_set1_epi32(value);
            v.x[1] = _mm_set1_epi32(value);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_set1_epi32(value);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = value;
        }
    #endif
    #endif
    
    return *this;
}

/** Destructor */
SIRE_ALWAYS_INLINE
MultiInt::~MultiInt()
{}

/** Comparison operator. This will return a MultiFloat with elements
    set to zero for each float that is not equal */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::compareEqual(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_cmpeq_epi32(v.x, other.v.x) );
        #else
            return MultiInt( _mm_cmpeq_epi32(v.x[0], other.v.x[0]),
                             _mm_cmpeq_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_cmpeq_epi32(v.x, other.v.x) );
    #else
        MultiInt ret;

        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] == other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Not equals comparison operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::compareNotEqual(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Less than comparison operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::compareLess(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            MultiInt ret;

            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIINT_BINONE : 0x0;
            }
            return ret;
        #else
            return MultiInt( _mm_cmplt_epi32(v.x[0], other.v.x[0]),
                             _mm_cmplt_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_cmplt_epi32(v.x, other.v.x) );
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Reintepret cast this MultiInt to a MultiFloat. This is only useful if you
    are going to use this to perform bitwise comparisons */
SIRE_ALWAYS_INLINE
MultiFloat MultiInt::reinterpretCastToFloat() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiFloat( *(reinterpret_cast<const __m256*>(&v.x)) );
        #else
            //allowable as v.x is 32bit aligned
            return MultiFloat( *(reinterpret_cast<const __m256*>(&v.x[0])) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( *(reinterpret_cast<const __m128*>(&v.x)) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = *(reinterpret_cast<const float*>(&v.a[i]));
        }
    
        return ret;
    #endif
    #endif
}

/** Greater than comparison operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::compareGreater(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_cmpgt_epi32(v.x, other.v.x) );
        #else
            return MultiInt( _mm_cmpgt_epi32(v.x[0], other.v.x[0]),
                             _mm_cmpgt_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_cmpgt_epi32(v.x, other.v.x) );
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Less than or equal comparison */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::compareLessEqual(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Greater than or equal comparison */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::compareGreaterEqual(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Return the number of values in the vector */
SIRE_ALWAYS_INLINE
int MultiInt::size()
{
    return MULTIFLOAT_SIZE;
}

/** Return the number of values in the vector */
SIRE_ALWAYS_INLINE
int MultiInt::count()
{
    return MULTIFLOAT_SIZE;
}

/** Addition operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator+(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_add_epi32(v.x, other.v.x) );
        #else
            return MultiInt( _mm_add_epi32(v.x[0], other.v.x[0]),
                             _mm_add_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_add_epi32(v.x, other.v.x) );
    #else
        MultiInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] + other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** Subtraction operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator-(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_sub_epi32(v.x, other.v.x) );
        #else
            return MultiInt( _mm_sub_epi32(v.x[0], other.v.x[0]),
                             _mm_sub_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_sub_epi32(v.x, other.v.x) );
    #else
        MultiInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] - other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** Multiplication operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator*(const MultiInt &other) const
{
    MultiInt ret;
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        ret.v.a[i] = v.a[i] * other.v.a[i];
    }
    return ret;
}

/** In-place addition operator */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator+=(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_add_epi32(v.x, other.v.x);
        #else
            v.x[0] = _mm_add_epi32(v.x[0], other.v.x[0]);
            v.x[1] = _mm_add_epi32(v.x[1], other.v.x[1]);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_add_epi32(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] += other.v.a[i];
        }
    #endif
    #endif

    return *this;
}

/** In-place subtraction operator */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator-=(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_sub_epi32(v.x, other.v.x);
        #else
            v.x[0] = _mm_sub_epi32(v.x[0], other.v.x[0]);
            v.x[1] = _mm_sub_epi32(v.x[1], other.v.x[1]);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_sub_epi32(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] -= other.v.a[i];
        }
    #endif
    #endif

    return *this;
}

/** In-place multiplication operator */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator*=(const MultiInt &other)
{
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        v.a[i] *= other.v.a[i];
    }

    return *this;
}

/** Bitwise logical "and" comparison */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::logicalAnd(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_and_si256(v.x, other.v.x) );
        #else
            return MultiInt( _mm_and_si128(v.x[0], other.v.x[0]),
                             _mm_and_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_and_si128(v.x, other.v.x) );
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(ret.v.a)[i] =
                reinterpret_cast<const unsigned int*>(other.v.a)[i] &
                reinterpret_cast<const unsigned int*>(v.a)[i];
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical "and" comparison */
SIRE_ALWAYS_INLINE
MultiFloat MultiFloat::logicalAnd(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiFloat(_mm256_and_ps(v.x,
                                *(reinterpret_cast<const __m256*>(&(other.v.x)))));
        #else
            //allowable as other.v.x is 32bit aligned
            return MultiFloat(_mm256_and_ps(v.x,
                                *(reinterpret_cast<const __m256*>(&(other.v.x[0])))));
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_and_ps(v.x, *(reinterpret_cast<const __m128*>(&(other.v.x)))) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(ret.v.a)[i] =
                reinterpret_cast<const unsigned int*>(other.v.a)[i] &
                reinterpret_cast<const unsigned int*>(v.a)[i];
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical "and not" */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::logicalAndNot(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_andnot_si256(v.x, other.v.x) );
        #else
            return MultiInt( _mm_andnot_si128(v.x[0], other.v.x[0]),
                             _mm_andnot_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_andnot_si128(v.x, other.v.x) );
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(ret.v.a)[i] =
                (~reinterpret_cast<const unsigned int*>(other.v.a)[i]) &
                reinterpret_cast<const unsigned int*>(v.a)[i];
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical "and not" (this is *this and (not other)) */
SIRE_ALWAYS_INLINE
MultiFloat MultiFloat::logicalAndNot(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            const __m256 val = *(reinterpret_cast<const __m256*>(&other.v.x));
            return MultiFloat( _mm256_andnot_ps(val, v.x) );
        #else
            //possible as other.v.x is 32bit aligned
            const __m256 val = *(reinterpret_cast<const __m256*>(&other.v.x[0]));
            return MultiFloat( _mm256_andnot_ps(val, v.x) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        const __m128 val = *(reinterpret_cast<const __m128*>(&other.v.x));
        return MultiFloat( _mm_andnot_ps(val, v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(ret.v.a)[i] =
                (~reinterpret_cast<const unsigned int*>(other.v.a)[i]) &
                reinterpret_cast<const unsigned int*>(v.a)[i];
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical or operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::logicalOr(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_or_si256(v.x, other.v.x) );
        #else
            return MultiInt( _mm_or_si128(v.x[0], other.v.x[0]),
                             _mm_or_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_or_si128(v.x, other.v.x) );
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(ret.v.a)[i] =
                reinterpret_cast<const unsigned int*>(other.v.a)[i] |
                reinterpret_cast<const unsigned int*>(v.a)[i];
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical xor */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::logicalXor(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiInt( _mm256_xor_si256(v.x, other.v.x) );
        #else
            return MultiInt( _mm_xor_si128(v.x[0], other.v.x[0]),
                             _mm_xor_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiInt( _mm_xor_si128(v.x, other.v.x) );
    #else
        MultiInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(ret.v.a)[i] =
                reinterpret_cast<const unsigned int*>(other.v.a)[i] ^
                reinterpret_cast<const unsigned int*>(v.a)[i];
        }
    
        return ret;
    #endif
    #endif
}

/** Logical not operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::logicalNot() const
{
    MultiInt ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        reinterpret_cast<unsigned int*>(ret.v.a)[i] =
            ~reinterpret_cast<const unsigned int*>(v.a)[i];
    }

    return ret;
}

/** Logical not operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator!() const
{
    return this->logicalNot();
}

/** Logical and operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator&(const MultiInt &other) const
{
    return this->logicalAnd(other);
}

/** Logical or operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator|(const MultiInt &other) const
{
    return this->logicalOr(other);
}

/** Logical xor operator */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::operator^(const MultiInt &other) const
{
    return this->logicalXor(other);
}

/** In place logical and */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator&=(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_and_si256(v.x, other.v.x);
        #else
            v.x[0] = _mm_and_si128(v.x[0], other.v.x[0]);
            v.x[1] = _mm_and_si128(v.x[1], other.v.x[1]);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_and_si128(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(v.a)[i] &=
                reinterpret_cast<const unsigned int*>(other.v.a)[i];
        }
    #endif
    #endif

    return *this;
}

/** In place logical and */
SIRE_ALWAYS_INLINE
MultiFloat& MultiFloat::operator&=(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_and_ps( v.x, *(reinterpret_cast<const __m256*>(&(other.v.x))) );
        #else
            //possible as other.v.x[0] is 32bit aligned
            v.x = _mm256_and_ps( v.x, *(reinterpret_cast<const __m256*>(&(other.v.x[0]))) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_and_ps( v.x, *(reinterpret_cast<const __m128*>(&(other.v.x))) );
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(v.a)[i] &=
                reinterpret_cast<const unsigned int*>(other.v.a)[i];
        }
    #endif
    #endif

    return *this;
}

/** In-place logical or */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator|=(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_or_si256(v.x, other.v.x);
        #else
            v.x[0] = _mm_or_si128(v.x[0], other.v.x[0]);
            v.x[1] = _mm_or_si128(v.x[1], other.v.x[1]);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_or_si128(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(v.a)[i] |=
                reinterpret_cast<const unsigned int*>(other.v.a)[i];
        }
    #endif
    #endif

    return *this;
}

/** In-place logical xor */
SIRE_ALWAYS_INLINE
MultiInt& MultiInt::operator^=(const MultiInt &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            v.x = _mm256_xor_si256(v.x, other.v.x);
        #else
            v.x[0] = _mm_xor_si128(v.x[0], other.v.x[0]);
            v.x[1] = _mm_xor_si128(v.x[1], other.v.x[1]);
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_xor_si128(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            reinterpret_cast<unsigned int*>(v.a)[i] ^=
                reinterpret_cast<const unsigned int*>(other.v.a)[i];
        }
    #endif
    #endif

    return *this;
}

/** Return the maximum vector between this and other */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::max(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
        }
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        #ifdef MULTIFLOAT_SSE4_IS_AVAILABLE
            return MultiInt( _mm_max_epi32(v.x, other.v.x) );
        #else
            MultiInt ret;
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
            }
            return ret;
        #endif
    #else
        MultiInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
        }
        return ret;
    #endif
    #endif
}

/** Return the minimum vector between this and other */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::min(const MultiInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
        }
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        #ifdef MULTIFLOAT_SSE4_IS_AVAILABLE
            return MultiInt( _mm_min_epi32(v.x, other.v.x) );
        #else
            MultiInt ret;
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
            }
            return ret;
        #endif
    #else
        MultiInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
        }
        return ret;
    #endif
    #endif
}

/** Rotate this vector. This moves each element one space to the left, moving the
    first element to the last element */
SIRE_ALWAYS_INLINE
MultiInt MultiInt::rotate() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiInt ret;
        
        for (int i=1; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i-1] = v.a[i];
        }
        
        ret.v.a[MULTIFLOAT_SIZE-1] = v.a[0];

        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        // there must be an SSE intrinsic to rotate left...
        return MultiInt( _mm_shuffle_epi32(v.x, _MM_SHUFFLE(0,3,2,1)) );
    #else
        MultiInt ret;
        
        for (int i=1; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i-1] = v.a[i];
        }
        
        ret.v.a[MULTIFLOAT_SIZE-1] = v.a[0];

        return ret;
    #endif
    #endif
}

/** Return the sum of all elements of this vector */
SIRE_ALWAYS_INLINE
qint32 MultiInt::sum() const
{
    qint32 sum = 0;
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        sum += v.a[i];
    }
    return sum;
}

/** Return the sum of all elements of this vector, using doubles for the sum */
SIRE_ALWAYS_INLINE
qint64 MultiInt::doubleSum() const
{
    qint64 sum = 0;
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        sum += v.a[i];
    }
    return sum;
}

#endif // #ifndef SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_TYPEINFO( SireMaths::MultiInt, Q_PRIMITIVE_TYPE );

SIRE_EXPOSE_CLASS( SireMaths::MultiInt )

SIRE_END_HEADER

#endif
