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

#ifndef SIREMATHS_MULTIFLOAT_H
#define SIREMATHS_MULTIFLOAT_H

#include "sireglobal.h"
#include <cmath>

#include <QVector>

#ifdef SIRE_HAS_CPP_11
    #include <functional>
#endif

#if defined(_MSC_VER)
#define _ALIGNED(x) __declspec(align(x))
#else
#if defined(__clang__)
#define _ALIGNED(x) __attribute__ ((aligned(x)))
#else
#if defined(__GNUC__)
#define _ALIGNED(x) __attribute__ ((aligned(x)))
#else
#define _ALIGNED(x)
#endif
#endif
#endif

//#undef SIRE_USE_AVX
//#define SIRE_USE_SSE 1
//#undef SIRE_USE_SSE

//#define MULTIFLOAT_CHECK_ALIGNMENT 1
//#undef MULTIFLOAT_CHECK_ALIGNMENT

#ifdef SIRE_USE_AVX
    #ifdef __AVX__
        #include <immintrin.h>   // CONDITIONAL_INCLUDE
        #define MULTIFLOAT_AVX_IS_AVAILABLE 1
        #undef MULTIFLOAT_SSE_IS_AVAILABLE

        #undef MULTIFLOAT_AVX2_IS_AVAILABLE

        #ifdef SIRE_USE_AVX2
            #ifdef __AVX2__
                #define MULTIFLOAT_AVX2_IS_AVAILABLE 1
            #endif
        #endif
    #else
    #ifdef __SSE2__
        #ifdef SIRE_USE_SSE4
           #ifdef __SSE4_1__
               #include <smmintrin.h>
               #define MULTIFLOAT_SSE4_IS_AVAILABLE 1
           #endif
        #endif
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
        #define MULTIFLOAT_SSE_IS_AVAILABLE 1
        #undef MULTIFLOAT_AVX_IS_AVAILABLE
    #else
        #undef SIRE_USE_SSE
        #undef MULTIFLOAT_SSE_IS_AVAILABLE
        #undef MULTIFLOAT_AVX_IS_AVAILABLE
    #endif
    #endif
#else
#ifdef SIRE_USE_SSE
    #ifdef __SSE2__
        #ifdef SIRE_USE_SSE4
           #ifdef __SSE4_1__
               #include <smmintrin.h>
               #define MULTIFLOAT_SSE4_IS_AVAILABLE 1
           #endif
        #endif
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
        #define MULTIFLOAT_SSE_IS_AVAILABLE 1
        #undef MULTIFLOAT_AVX_IS_AVAILABLE
    #else
        #undef SIRE_USE_SSE
        #undef MULTIFLOAT_SSE_IS_AVAILABLE
        #undef MULTIFLOAT_AVX_IS_AVAILABLE
    #endif
#else
    #undef MULTIFLOAT_SSE_IS_AVAILABLE
    #undef MULTIFLOAT_AVX_IS_AVAILABLE
#endif
#endif

SIRE_BEGIN_HEADER

namespace SireMaths
{

class MultiFloat;
class MultiFixed;
class MultiDouble;
class MultiUInt;
class MultiInt;

MultiFloat cos(const MultiFloat &val);
MultiFloat sin(const MultiFloat &val);
MultiFloat log(const MultiFloat &val);
MultiFloat exp(const MultiFloat &val);

void sincos(const MultiFloat &val, MultiFloat &sinval, MultiFloat &cosval);

/** This class provides a vectorised float. This represents
    a single vector of floats on the compiled machine, e.g.
    4 floats if we use SSE2, 8 floats for AVX
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT MultiFloat
{

friend MultiFloat cos(const MultiFloat &val);
friend MultiFloat sin(const MultiFloat &val);
friend MultiFloat log(const MultiFloat &val);
friend MultiFloat exp(const MultiFloat &val);

friend void sincos(const MultiFloat &val, MultiFloat &sinval, MultiFloat &cosval);

public:
    MultiFloat();
    
    MultiFloat(float value);
    
    MultiFloat(const float *array, int size);
    MultiFloat(const float *array, const MultiInt &indicies);
    
    MultiFloat(const QVector<float> &array);
    MultiFloat(const QVector<double> &array);
    
    MultiFloat(const MultiDouble &other);
    MultiFloat(const MultiFloat &other);
    MultiFloat(const MultiInt &other);
    
    #ifdef SIRE_HAS_CPP_11
        MultiFloat(const std::function<float ()> &func);
    #endif
    
    ~MultiFloat();
    
    bool isAligned() const;
    
    static QVector<MultiFloat> fromArray(const QVector<double> &array);
    static QVector<MultiFloat> fromArray(const QVector<float> &array);
    
    static QVector<MultiFloat> fromArray(const double *array, int size);
    static QVector<MultiFloat> fromArray(const float *array, int size);
    
    static QVector<float> toArray(const QVector<MultiFloat> &array);
    static QVector<double> toDoubleArray(const QVector<MultiFloat> &array);
    
    MultiFloat& operator=(const MultiFloat &other);
    MultiFloat& operator=(const MultiDouble &other);
    MultiFloat& operator=(const MultiInt &other);
    MultiFloat& operator=(float value);
    
    bool operator==(const MultiFloat &other) const;
    bool operator!=(const MultiFloat &other) const;
    
    bool operator<(const MultiFloat &other) const;
    bool operator>(const MultiFloat &other) const;
    
    bool operator<=(const MultiFloat &other) const;
    bool operator>=(const MultiFloat &other) const;
    
    MultiFloat compareEqual(const MultiFloat &other) const;
    MultiFloat compareNotEqual(const MultiFloat &other) const;

    MultiFloat compareLess(const MultiFloat &other) const;
    MultiFloat compareGreater(const MultiFloat &other) const;
    
    MultiFloat compareLessEqual(const MultiFloat &other) const;
    MultiFloat compareGreaterEqual(const MultiFloat &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    QString toBinaryString() const;
    
    static int size();
    static int count();
    
    float operator[](int i) const;
    
    float at(int i) const;
    float getitem(int i) const;
    
    void set(int i, float value);
    float get(int i) const;
    
    void quickSet(int i, float value);
    
    MultiFloat operator-() const;
    
    MultiFloat operator+(const MultiFloat &other) const;
    MultiFloat operator-(const MultiFloat &other) const;
    MultiFloat operator*(const MultiFloat &other) const;
    MultiFloat operator/(const MultiFloat &other) const;
    
    MultiFloat& operator+=(const MultiFloat &other);
    MultiFloat& operator-=(const MultiFloat &other);
    MultiFloat& operator*=(const MultiFloat &other);
    MultiFloat& operator/=(const MultiFloat &other);
    
    MultiFloat operator!() const;
    MultiFloat operator&(const MultiFloat &other) const;
    MultiFloat operator|(const MultiFloat &other) const;
    MultiFloat operator^(const MultiFloat &other) const;

    MultiFloat& operator&=(const MultiFloat &other);
    MultiFloat& operator|=(const MultiFloat &other);
    MultiFloat& operator^=(const MultiFloat &other);

    MultiFloat& operator&=(const MultiInt &other);
    MultiFloat& operator&=(const MultiUInt &other);

    MultiFloat logicalNot() const;
    
    MultiFloat logicalAnd(const MultiFloat &other) const;
    MultiFloat logicalAndNot(const MultiFloat &other) const;
    
    MultiFloat logicalOr(const MultiFloat &other) const;
    MultiFloat logicalXor(const MultiFloat &other) const;
    
    MultiFloat logicalAnd(const MultiUInt &other) const;
    MultiFloat logicalAnd(const MultiInt &other) const;
    
    MultiFloat logicalAndNot(const MultiInt &other) const;
    MultiFloat logicalAndNot(const MultiUInt &other) const;
    
    MultiFloat& multiplyAdd(const MultiFloat &val0, const MultiFloat &val1);
    
    MultiFloat max(const MultiFloat &other) const;
    MultiFloat min(const MultiFloat &other) const;
    
    float max() const;
    float min() const;
    
    MultiFloat reciprocal() const;
    MultiFloat reciprocal_approx() const;
    MultiFloat reciprocal_approx_nr() const;
    
    MultiFloat sqrt() const;
    MultiFloat sqrt_approx() const;
    MultiFloat sqrt_approx_nr() const;
    
    MultiFloat rsqrt() const;
    MultiFloat rsqrt_approx() const;
    MultiFloat rsqrt_approx_nr() const;
    
    MultiFloat rotate() const;
    
    MultiFloat abs() const;
    
    bool isBinaryZero() const;
    bool isNotBinaryZero() const;
    bool hasBinaryZero() const;
    
    bool isBinaryOne() const;
    bool isNotBinaryOne() const;
    bool hasBinaryOne() const;
    
    float sum() const;
    double doubleSum() const;

    static void swap(MultiFloat &f0, int idx0, MultiFloat &f1, int idx1);

private:
    /* Make the other Multi?? classes friends, so that they
       can read the vector data directly */
    friend class MultiDouble;
    friend class MultiFixed;
    friend class MultiInt;
    friend class MultiUInt;

    static void assertAligned(const void *ptr, size_t size);

    #ifndef SIRE_SKIP_INLINE_FUNCTIONS
        #ifndef MULTIFLOAT_CHECK_ALIGNMENT
            void assertAligned(){}
        #endif

        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            _ALIGNED(32) union
            {
                __m256 x;
                float a[8];
            } v;
            #define MULTIFLOAT_SIZE 8
        
            MultiFloat(__m256 avx_val)
            {
                v.x = avx_val;
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
                __m128 x;
                float a[4];
            } v;
            #define MULTIFLOAT_SIZE 4

            #define MULTIFLOAT_BINONE getBinaryOne()

            MultiFloat(__m128 sse_val)
            {
                v.x = sse_val;
            }

            static float getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return *(reinterpret_cast<const float*>(&x));
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
                float a[8];
            } v;
            #define MULTIFLOAT_SIZE 8
            #define MULTIFLOAT_BINONE getBinaryOne()

            static float getBinaryOne()
            {
                const quint32 x = 0xFFFFFFFF;
                return *(reinterpret_cast<const float*>(&x));
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
    #else
        #define MULTIFLOAT_SIZE 16
    #endif //#ifndef SIRE_SKIP_INLINE_FUNCTIONS

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

static const MultiFloat MULTIFLOAT_ONE(1);
static const MultiFloat MULTIFLOAT_NEGATIVE_ONE(-1);

static quint32 MULTIFLOAT_POS_MASK_INT = 0x7FFFFFFF;
static const MultiFloat MULTIFLOAT_POS_MASK(
                            *(reinterpret_cast<const float*>(&MULTIFLOAT_POS_MASK_INT)) );

/** Constructor. This creates a MultiFloat with an undefined initial state */
inline
MultiFloat::MultiFloat()
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_set1_ps(0);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_set1_ps(0);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = 0;
        }
    #endif
    #endif
}

/** Construct a MultiFloat with all values equal to 'val' */
inline
MultiFloat::MultiFloat(float val)
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_set1_ps(val);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_set1_ps(val);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = val;
        }
    #endif
    #endif
}

#ifdef SIRE_HAS_CPP_11
    /** Construct from the passed function that generates floats */
    inline MultiFloat::MultiFloat(const std::function<float ()> &generator)
    {
        assertAligned();
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            v.x = _mm256_set_ps(generator(), generator(), generator(), generator(),
                                generator(), generator(), generator(), generator());
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            v.x = _mm_set_ps(generator(), generator(), generator(), generator());
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = generator();
            }
        #endif
        #endif
    }
#endif

/** Copy constructor */
inline
MultiFloat::MultiFloat(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = other.v.x;
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

/** Assignment operator */
inline
MultiFloat& MultiFloat::operator=(const MultiFloat &other)
{
    if (this != &other)
    {
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            v.x = other.v.x;
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
inline
MultiFloat& MultiFloat::operator=(float value)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_set1_ps(value);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_set1_ps(value);
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
inline
MultiFloat::~MultiFloat()
{}

/** Return the ith value in the MultiFloat. This is a fast function
    that performs no bounds checking! */
inline float MultiFloat::operator[](int i) const
{
    return v.a[i];
}

/** Comparison operator. This will return a MultiFloat with elements
    set to zero for each float that is not equal */
inline
MultiFloat MultiFloat::compareEqual(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_cmp_ps(v.x, other.v.x, _CMP_EQ_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_cmpeq_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;

        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] == other.v.a[i]) ? MULTIFLOAT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Not equals comparison operator */
inline
MultiFloat MultiFloat::compareNotEqual(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_cmp_ps(v.x, other.v.x, _CMP_NEQ_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_cmpneq_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIFLOAT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Less than comparison operator */
inline
MultiFloat MultiFloat::compareLess(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_cmp_ps(v.x, other.v.x, _CMP_LT_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_cmplt_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIFLOAT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Greater than comparison operator */
inline
MultiFloat MultiFloat::compareGreater(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_cmp_ps(v.x, other.v.x, _CMP_GT_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_cmpgt_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIFLOAT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Less than or equal comparison */
inline
MultiFloat MultiFloat::compareLessEqual(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_cmp_ps(v.x, other.v.x, _CMP_LE_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_cmple_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIFLOAT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Greater than or equal comparison */
inline
MultiFloat MultiFloat::compareGreaterEqual(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_cmp_ps(v.x, other.v.x, _CMP_GE_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_cmpge_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIFLOAT_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
}

/** Return the number of values in the vector */
inline
int MultiFloat::size()
{
    return MULTIFLOAT_SIZE;
}

/** Return the number of values in the vector */
inline
int MultiFloat::count()
{
    return MULTIFLOAT_SIZE;
}

/** Addition operator */
inline
MultiFloat MultiFloat::operator+(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_add_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_add_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] + other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** Subtraction operator */
inline
MultiFloat MultiFloat::operator-(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_sub_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_sub_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] - other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** Multiplication operator */
inline
MultiFloat MultiFloat::operator*(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_mul_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_mul_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] * other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** Division operator */
inline
MultiFloat MultiFloat::operator/(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_div_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_div_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] / other.v.a[i];
        }
        return ret;
    #endif
    #endif
}

/** In-place addition operator */
inline
MultiFloat& MultiFloat::operator+=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_add_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_add_ps(v.x, other.v.x);
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
inline
MultiFloat& MultiFloat::operator-=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_sub_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_sub_ps(v.x, other.v.x);
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
inline
MultiFloat& MultiFloat::operator*=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_mul_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_mul_ps(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] *= other.v.a[i];
        }
    #endif
    #endif

    return *this;
}

/** In-place division operator */
inline
MultiFloat& MultiFloat::operator/=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_div_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_div_ps(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] /= other.v.a[i];
        }
    #endif
    #endif

    return *this;
}

/** Negative operator */
inline MultiFloat MultiFloat::operator-() const
{
    return this->operator*(MULTIFLOAT_NEGATIVE_ONE);
}

/** Bitwise logical "and" comparison */
inline
MultiFloat MultiFloat::logicalAnd(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_and_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_and_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                ret_char_v[j] = char_v[j] & other_char_v[j];
            }
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical "and not" (this is *this and (not other)) */
inline
MultiFloat MultiFloat::logicalAndNot(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_andnot_ps(other.v.x, v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_andnot_ps(other.v.x, v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                ret_char_v[j] = !(char_v[j] & other_char_v[j]);
            }
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical or operator */
inline
MultiFloat MultiFloat::logicalOr(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_or_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_or_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                ret_char_v[j] = char_v[j] | other_char_v[j];
            }
        }
    
        return ret;
    #endif
    #endif
}

/** Bitwise logical xor */
inline
MultiFloat MultiFloat::logicalXor(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_xor_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_xor_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                ret_char_v[j] = char_v[j] ^ other_char_v[j];
            }
        }
    
        return ret;
    #endif
    #endif
}

/** Logical not operator */
inline
MultiFloat MultiFloat::logicalNot() const
{
    MultiFloat ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));

        for (unsigned int j=0; j<sizeof(float); ++j)
        {
            ret_char_v[j] = !char_v[j];
        }
    }

    return ret;
}

/** Logical not operator */
inline
MultiFloat MultiFloat::operator!() const
{
    return this->logicalNot();
}

/** Logical and operator */
inline
MultiFloat MultiFloat::operator&(const MultiFloat &other) const
{
    return this->logicalAnd(other);
}

/** Logical or operator */
inline
MultiFloat MultiFloat::operator|(const MultiFloat &other) const
{
    return this->logicalOr(other);
}

/** Logical xor operator */
inline
MultiFloat MultiFloat::operator^(const MultiFloat &other) const
{
    return this->logicalXor(other);
}

/** In place logical and */
inline
MultiFloat& MultiFloat::operator&=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_and_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_and_ps(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                char_v[j] &= other_char_v[j];
            }
        }
    #endif
    #endif

    return *this;
}

/** In-place logical or */
inline
MultiFloat& MultiFloat::operator|=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_or_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_or_ps(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                char_v[j] |= other_char_v[j];
            }
        }
    #endif
    #endif

    return *this;
}

/** In-place logical xor */
inline
MultiFloat& MultiFloat::operator^=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_xor_ps(v.x, other.v.x);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_xor_ps(v.x, other.v.x);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(float); ++j)
            {
                char_v[j] ^= other_char_v[j];
            }
        }
    #endif
    #endif

    return *this;
}

/** Multiply val0 and val1 and add it onto this value */
inline
MultiFloat& MultiFloat::multiplyAdd(const MultiFloat &v0, const MultiFloat &v1)
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x = _mm256_add_ps(v.x, _mm256_mul_ps(v0.v.x, v1.v.x));
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_add_ps(v.x, _mm_mul_ps(v0.v.x, v1.v.x));
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] += v0.v.a[i] * v1.v.a[i];
        }
    #endif
    #endif

    return *this;
}

/** Return the maximum vector between this and other */
inline
MultiFloat MultiFloat::max(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_max_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_max_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
        }
        return ret;
    #endif
    #endif
}

/** Return the minimum vector between this and other */
inline
MultiFloat MultiFloat::min(const MultiFloat &other) const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_min_ps(v.x, other.v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_min_ps(v.x, other.v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
        }
        return ret;
    #endif
    #endif
}

/** Return the reciprocal of this vector */
inline
MultiFloat MultiFloat::reciprocal() const
{
    return MULTIFLOAT_ONE.operator/(*this);
}

/** Return a poor approximation of the reciprocal of the vector (about 12 bits of precision) */
inline
MultiFloat MultiFloat::reciprocal_approx() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_rcp_ps(v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_rcp_ps(v.x) );
    #else
        return this->reciprocal();
    #endif
    #endif
}

/** Return a good approximation of the reciprocal of the vector (the poor approximation
    refined using one step of Newton Raphson) */
inline
MultiFloat MultiFloat::reciprocal_approx_nr() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        //get the approximation
        __m256 a = _mm256_rcp_ps(v.x);

        //now use one step of NR to refine the result
        // 1/x = a[ 2 - a x ] where a is the approximation
        __m256 tmp = _mm256_mul_ps(a, v.x);
        __m256 two = _mm256_set1_ps(2.0);
        tmp = _mm256_sub_ps(two, tmp);
        a = _mm256_mul_ps(a, tmp);
    
        return MultiFloat(a);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        //get the approximation
        __m128 a = _mm_rcp_ps(v.x);

        //now use one step of NR to refine the result
        // 1/x = a[ 2 - a x ] where a is the approximation
        __m128 tmp = _mm_mul_ps(a, v.x);
        __m128 two = _mm_set1_ps(2.0);
        tmp = _mm_sub_ps(two, tmp);
        a = _mm_mul_ps(a, tmp);
    
        return MultiFloat(a);
    #else
        return this->reciprocal();
    #endif
    #endif
}

/** Return the square root of this vector */
inline
MultiFloat MultiFloat::sqrt() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_sqrt_ps(v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_sqrt_ps(v.x) );
    #else
        MultiFloat ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::sqrt(v.a[i]);
        }
        return ret;
    #endif
    #endif
}

/** Return the approximate square root, highly approximated but very fast.
    Has about 12 bits of precision in the mantissa (is x * rsqrt(x)) */
inline
MultiFloat MultiFloat::sqrt_approx() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        __m256 r_sqrt = _mm256_rsqrt_ps(v.x);
        return MultiFloat( _mm256_mul_ps( v.x, r_sqrt ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        __m128 r_sqrt = _mm_rsqrt_ps(v.x);
        return MultiFloat( _mm_mul_ps( v.x, r_sqrt ) );
    #else
        return this->sqrt();
    #endif
    #endif
}

/** Return a good approximation of the square root (this is the poor approximation
    refined using a single step of Newton Raphson) */
inline
MultiFloat MultiFloat::sqrt_approx_nr() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        //calculate sqrt(x) as x * 1 / sqrt(x)
        __m256 a = _mm256_rsqrt_ps(v.x);
        a = _mm256_mul_ps(v.x, a);

        //now use one step of NR to refine the result
        // sqrt(x) = a - [ (a^2 - x) / 2a ] where a is the approximation
        __m256 tmp = _mm256_mul_ps(a, a);
        tmp = _mm256_sub_ps(tmp, v.x);
        __m256 two_a = _mm256_add_ps(a, a);
        tmp = _mm256_div_ps(tmp, two_a);
        a = _mm256_sub_ps(a, tmp);
    
        return MultiFloat(a);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        //calculate sqrt(x) as x * 1 / sqrt(x)
        __m128 a = _mm_rsqrt_ps(v.x);
        a = _mm_mul_ps(v.x, a);

        //now use one step of NR to refine the result
        // sqrt(x) = a - [ (a^2 - x) / 2a ] where a is the approximation
        __m128 tmp = _mm_mul_ps(a, a);
        tmp = _mm_sub_ps(tmp, v.x);
        __m128 two_a = _mm_add_ps(a, a);
        tmp = _mm_div_ps(tmp, two_a);
        a = _mm_sub_ps(a, tmp);
    
        return MultiFloat(a);
    #else
        return this->sqrt();
    #endif
    #endif
}

/** Return the recipical square root of this vector */
inline
MultiFloat MultiFloat::rsqrt() const
{
    return MULTIFLOAT_ONE.operator/(this->sqrt());
}

/** Return an approximation of the reciprical square root of this vector (only good
    for about 12 bit of the mantissa) */
inline
MultiFloat MultiFloat::rsqrt_approx() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiFloat( _mm256_rsqrt_ps(v.x) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiFloat( _mm_rsqrt_ps(v.x) );
    #else
        return MULTIFLOAT_ONE.operator/(this->sqrt());
    #endif
    #endif
}

/** Return a good approximation of the reciprical square root (this poor approximation
    refined using a single Newton Raphson step) */
inline
MultiFloat MultiFloat::rsqrt_approx_nr() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        //get the approximation
        __m256 a = _mm256_rsqrt_ps(v.x);

        //now use one step of NR to refine the result
        // 1/x = 0.5 a[ 3 - x a^2 ] where a is the approximation
        __m256 tmp = _mm256_mul_ps(a, v.x);
        tmp = _mm256_mul_ps(a, tmp);
        __m256 three = _mm256_set1_ps(3.0);
        tmp = _mm256_sub_ps(three, tmp);
        a = _mm256_mul_ps(a, tmp);
        __m256 half = _mm256_set1_ps(0.5);
        a = _mm256_mul_ps(a, half);
    
        return MultiFloat(a);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        //get the approximation
        __m128 a = _mm_rsqrt_ps(v.x);

        //now use one step of NR to refine the result
        // 1/x = 0.5 a[ 3 - x a^2 ] where a is the approximation
        __m128 tmp = _mm_mul_ps(a, v.x);
        tmp = _mm_mul_ps(a, tmp);
        __m128 three = _mm_set1_ps(3.0);
        tmp = _mm_sub_ps(three, tmp);
        a = _mm_mul_ps(a, tmp);
        __m128 half = _mm_set1_ps(0.5);
        a = _mm_mul_ps(a, half);
    
        return MultiFloat(a);
    #else
        return this->rsqrt();
    #endif
    #endif
}

/** Rotate this vector. This moves each element one space to the left, moving the
    first element to the last element */
inline
MultiFloat MultiFloat::rotate() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        __m256 tmp =  _mm256_permute_ps(v.x, _MM_SHUFFLE ( 0,3,2,1 ));
        return MultiFloat( _mm256_blend_ps(tmp, _mm256_permute2f128_ps ( tmp,tmp,1 ), 136) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        // there must be an SSE intrinsic to rotate left...
        return MultiFloat( _mm_shuffle_ps(v.x, v.x, _MM_SHUFFLE(0,3,2,1)) );
    #else
        MultiFloat ret;
        
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
inline
float MultiFloat::sum() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return v.a[0] + v.a[1] + v.a[2] + v.a[3] +
               v.a[4] + v.a[5] + v.a[6] + v.a[7];
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return v.a[0] + v.a[1] + v.a[2] + v.a[3];
    #else
        float sum = 0;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            sum += v.a[i];
        }
        return sum;
    #endif
    #endif
}

/** Return the sum of all elements of this vector, using doubles for the sum */
inline
double MultiFloat::doubleSum() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return double(v.a[0]) + double(v.a[1]) + double(v.a[2]) + double(v.a[3]) +
               double(v.a[4]) + double(v.a[5]) + double(v.a[6]) + double(v.a[7]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return double(v.a[0]) + double(v.a[1]) + double(v.a[2]) + double(v.a[3]);
    #else
        double sum = 0;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            sum += double(v.a[i]);
        }
        return sum;
    #endif
    #endif
}

/** Return the absolute (positive) value of the floats */
inline MultiFloat MultiFloat::abs() const
{
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        MultiFloat ret;
        ret.v.x = _mm256_and_ps(v.x, MULTIFLOAT_POS_MASK.v.x);
        return ret;
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        MultiFloat ret;
        ret.v.x = _mm_and_ps(v.x, MULTIFLOAT_POS_MASK.v.x);
        return ret;
    #else
        MultiFloat ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] < 0 ? -v.a[i] : v.a[i]);
        }
    
        return ret;
    #endif
    #endif
}

/** Quick function to set the value of the ith element */
inline void MultiFloat::quickSet(int i, float value)
{
    v.a[i] = value;
}

#endif // #ifndef SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_TYPEINFO( SireMaths::MultiFloat, Q_PRIMITIVE_TYPE );

SIRE_EXPOSE_CLASS( SireMaths::MultiFloat )

SIRE_END_HEADER

#endif
