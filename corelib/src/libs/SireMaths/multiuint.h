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

#ifndef SIREMATHS_MULTIUINT_H
#define SIREMATHS_MULTIUINT_H

#include "multifloat.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

/** This class provides a vectorised 32bit unsigned integer. This represents
    a single vector of integers on the compiled machine, e.g.
    4 integers if we use SSE2, 8 integers for AVX/AVX2
    (note that AVX represents it as two SSE vectors, while AVX2 
     uses a single AVX vector)
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT MultiUInt
{
public:
    MultiUInt();
    
    MultiUInt(quint32 value);
    
    MultiUInt(const quint32 *array, int size);
    MultiUInt(const QVector<quint32> &array);
    
    MultiUInt(const MultiUInt &other);
    
    ~MultiUInt();
    
    bool isAligned() const;
    
    static QVector<MultiUInt> fromArray(const QVector<quint32> &array);
    static QVector<MultiUInt> fromArray(const quint32 *array, int size);
    
    static QVector<quint32> toArray(const QVector<MultiUInt> &array);
    
    MultiUInt& operator=(const MultiUInt &other);
    MultiUInt& operator=(quint32 value);
    
    bool operator==(const MultiUInt &other) const;
    bool operator!=(const MultiUInt &other) const;
    
    bool operator<(const MultiUInt &other) const;
    bool operator>(const MultiUInt &other) const;
    
    bool operator<=(const MultiUInt &other) const;
    bool operator>=(const MultiUInt &other) const;
    
    MultiUInt compareEqual(const MultiUInt &other) const;
    MultiUInt compareNotEqual(const MultiUInt &other) const;

    MultiUInt compareLess(const MultiUInt &other) const;
    MultiUInt compareGreater(const MultiUInt &other) const;
    
    MultiUInt compareLessEqual(const MultiUInt &other) const;
    MultiUInt compareGreaterEqual(const MultiUInt &other) const;
    
    MultiFloat reinterpretCastToFloat() const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    QString toBinaryString() const;
    
    static int size();
    static int count();
    
    quint32 operator[](int i) const;
    
    void set(int i, quint32 value);
    quint32 get(int i) const;
    
    quint32 at(int i) const;
    quint32 getitem(int i) const;
    
    MultiUInt operator+(const MultiUInt &other) const;
    MultiUInt operator-(const MultiUInt &other) const;
    
    MultiUInt& operator+=(const MultiUInt &other);
    MultiUInt& operator-=(const MultiUInt &other);
    
    MultiUInt operator!() const;
    MultiUInt operator&(const MultiUInt &other) const;
    MultiUInt operator|(const MultiUInt &other) const;
    MultiUInt operator^(const MultiUInt &other) const;

    MultiUInt& operator&=(const MultiUInt &other);
    MultiUInt& operator|=(const MultiUInt &other);
    MultiUInt& operator^=(const MultiUInt &other);

    MultiUInt logicalNot() const;
    
    MultiUInt logicalAnd(const MultiUInt &other) const;
    MultiUInt logicalAndNot(const MultiUInt &other) const;
    
    MultiUInt logicalOr(const MultiUInt &other) const;
    MultiUInt logicalXor(const MultiUInt &other) const;
    
    MultiUInt max(const MultiUInt &other) const;
    MultiUInt min(const MultiUInt &other) const;
    
    MultiUInt rotate() const;
    
    bool isBinaryZero() const;
    bool isNotBinaryZero() const;
    bool hasBinaryZero() const;
    
    bool isBinaryOne() const;
    bool isNotBinaryOne() const;
    bool hasBinaryOne() const;
    
    quint32 sum() const;
    quint64 doubleSum() const;

private:
    /* Make the other Multi?? classes friends, so that they
       can read the vector data directly */
    friend class MultiFloat;
    friend class MultiInt;

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
                    quint32 a[8];
                } v;
    
                MultiUInt(__m256i val)
                {
                    v.x = val;
                }
            #else
                _ALIGNED(32) union
                {
                    __m128i x[2];
                    quint32 a[8];
                } v;
        
                MultiUInt(__m128i val0, __m128i val1)
                {
                    v.x[0] = val0;
                    v.x[1] = val1;
                }
            #endif

            #define MULTIUINT_BINONE getBinaryOne()

            static quint32 getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return x;
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
                quint32 a[4];
            } v;

            MultiUInt(__m128i sse_val)
            {
                v.x = sse_val;
            }

            #define MULTIUINT_BINONE getBinaryOne()

            static quint32 getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return x;
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
                quint32 a[MULTIFLOAT_SIZE];
            } v;
            #define MULTIUINT_BINONE getBinaryOne()

            static quint32 getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return x;
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
MultiUInt::MultiUInt()
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
MultiUInt::MultiUInt(quint32 uval)
{
    assertAligned();

    qint32 val = *(reinterpret_cast<qint32*>(&uval));

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
MultiUInt::MultiUInt(const MultiUInt &other)
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

/** Return the ith value in the MultiUInt - note that this is
    a quick function that does no bounds checking */
SIRE_ALWAYS_INLINE quint32 MultiUInt::operator[](int i) const
{
    return *(reinterpret_cast<const quint32*>(&(v.a[i])));
}

/** Assignment operator */
SIRE_ALWAYS_INLINE
MultiUInt& MultiUInt::operator=(const MultiUInt &other)
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
MultiUInt& MultiUInt::operator=(quint32 uvalue)
{
    qint32 value = *(reinterpret_cast<qint32*>(&uvalue));

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
MultiUInt::~MultiUInt()
{}

/** Comparison operator. This will return a MultiFloat with elements
    set to zero for each float that is not equal */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::compareEqual(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_cmpeq_epi32(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_cmpeq_epi32(v.x[0], other.v.x[0]),
                              _mm_cmpeq_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_cmpeq_epi32(v.x, other.v.x) );
    #else
        MultiUInt ret;

        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] == other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Not equals comparison operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::compareNotEqual(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Less than comparison operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::compareLess(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            MultiUInt ret;

            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
            }
            return ret;
        #else
            MultiUInt ret;

            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
            }
            return ret;
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiUInt ret;

        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
        return ret;
    #else
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Reintepret cast this MultiInt to a MultiFloat. This is only useful if you
    are going to use this to perform bitwise comparisons */
SIRE_ALWAYS_INLINE
MultiFloat MultiUInt::reinterpretCastToFloat() const
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
MultiUInt MultiUInt::compareGreater(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            MultiUInt ret;
        
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
            }
        
            return ret;
        #else
            MultiUInt ret;
        
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
            }
        
            return ret;
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Less than or equal comparison */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::compareLessEqual(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Greater than or equal comparison */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::compareGreaterEqual(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #else
        MultiUInt ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIUINT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Return the number of values in the vector */
SIRE_ALWAYS_INLINE
int MultiUInt::size()
{
    return MULTIFLOAT_SIZE;
}

/** Return the number of values in the vector */
SIRE_ALWAYS_INLINE
int MultiUInt::count()
{
    return MULTIFLOAT_SIZE;
}

/** Addition operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::operator+(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_add_epi32(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_add_epi32(v.x[0], other.v.x[0]),
                              _mm_add_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_add_epi32(v.x, other.v.x) );
    #else
        MultiUInt ret;
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
MultiUInt MultiUInt::operator-(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_sub_epi32(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_sub_epi32(v.x[0], other.v.x[0]),
                              _mm_sub_epi32(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_sub_epi32(v.x, other.v.x) );
    #else
        MultiUInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] - other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** In-place addition operator */
SIRE_ALWAYS_INLINE
MultiUInt& MultiUInt::operator+=(const MultiUInt &other)
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
MultiUInt& MultiUInt::operator-=(const MultiUInt &other)
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

/** Bitwise logical "and" comparison */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::logicalAnd(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_and_si256(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_and_si128(v.x[0], other.v.x[0]),
                              _mm_and_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_and_si128(v.x, other.v.x) );
    #else
        MultiUInt ret;
    
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
MultiFloat MultiFloat::logicalAnd(const MultiUInt &other) const
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
MultiUInt MultiUInt::logicalAndNot(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_andnot_si256(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_andnot_si128(v.x[0], other.v.x[0]),
                              _mm_andnot_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_andnot_si128(v.x, other.v.x) );
    #else
        MultiUInt ret;
    
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
MultiFloat MultiFloat::logicalAndNot(const MultiUInt &other) const
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
MultiUInt MultiUInt::logicalOr(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_or_si256(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_or_si128(v.x[0], other.v.x[0]),
                              _mm_or_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_or_si128(v.x, other.v.x) );
    #else
        MultiUInt ret;
    
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
MultiUInt MultiUInt::logicalXor(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        #ifdef MULTIFLOAT_AVX2_IS_AVAILABLE
            return MultiUInt( _mm256_xor_si256(v.x, other.v.x) );
        #else
            return MultiUInt( _mm_xor_si128(v.x[0], other.v.x[0]),
                              _mm_xor_si128(v.x[1], other.v.x[1]) );
        #endif
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiUInt( _mm_xor_si128(v.x, other.v.x) );
    #else
        MultiUInt ret;
    
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
MultiUInt MultiUInt::logicalNot() const
{
    MultiUInt ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        reinterpret_cast<unsigned int*>(ret.v.a)[i] =
            ~reinterpret_cast<const unsigned int*>(v.a)[i];
    }

    return ret;
}

/** Logical not operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::operator!() const
{
    return this->logicalNot();
}

/** Logical and operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::operator&(const MultiUInt &other) const
{
    return this->logicalAnd(other);
}

/** Logical or operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::operator|(const MultiUInt &other) const
{
    return this->logicalOr(other);
}

/** Logical xor operator */
SIRE_ALWAYS_INLINE
MultiUInt MultiUInt::operator^(const MultiUInt &other) const
{
    return this->logicalXor(other);
}

/** In place logical and */
SIRE_ALWAYS_INLINE
MultiUInt& MultiUInt::operator&=(const MultiUInt &other)
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
MultiFloat& MultiFloat::operator&=(const MultiUInt &other)
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
MultiUInt& MultiUInt::operator|=(const MultiUInt &other)
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
MultiUInt& MultiUInt::operator^=(const MultiUInt &other)
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
MultiUInt MultiUInt::max(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiUInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
        }
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        #ifdef MULTIFLOAT_SSE4_IS_AVAILABLE
            return MultiUInt( _mm_max_epi32(v.x, other.v.x) );
        #else
            MultiUInt ret;
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
            }
            return ret;
        #endif
    #else
        MultiUInt ret;
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
MultiUInt MultiUInt::min(const MultiUInt &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiUInt ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
        }
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        #ifdef MULTIFLOAT_SSE4_IS_AVAILABLE
            return MultiUInt( _mm_min_epi32(v.x, other.v.x) );
        #else
            MultiUInt ret;
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
            }
            return ret;
        #endif
    #else
        MultiUInt ret;
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
MultiUInt MultiUInt::rotate() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiUInt ret;
        
        for (int i=1; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i-1] = v.a[i];
        }
        
        ret.v.a[MULTIFLOAT_SIZE-1] = v.a[0];

        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        // there must be an SSE intrinsic to rotate left...
        return MultiUInt( _mm_shuffle_epi32(v.x, _MM_SHUFFLE(0,3,2,1)) );
    #else
        MultiUInt ret;
        
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
quint32 MultiUInt::sum() const
{
    quint32 sum = 0;
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        sum += v.a[i];
    }
    return sum;
}

/** Return the sum of all elements of this vector, using doubles for the sum */
SIRE_ALWAYS_INLINE
quint64 MultiUInt::doubleSum() const
{
    quint64 sum = 0;
    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        sum += v.a[i];
    }
    return sum;
}

#endif // #ifndef SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireMaths::MultiUInt )

SIRE_END_HEADER

#endif
