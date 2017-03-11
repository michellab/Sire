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

#ifndef SIREMATHS_MULTIDOUBLE_H
#define SIREMATHS_MULTIDOUBLE_H

#include "SireMaths/multifloat.h"

#ifdef SIRE_HAS_CPP_11
    #include <functional>
#endif

SIRE_BEGIN_HEADER

namespace SireMaths
{

class MultiDouble;

MultiDouble cos(const MultiDouble &val);
MultiDouble sin(const MultiDouble &val);
MultiDouble log(const MultiDouble &val);
MultiDouble exp(const MultiDouble &val);

void sincos(const MultiDouble &val, MultiDouble &sinval, MultiDouble &cosval);

/** This class provides a vectorised double. This represents
    two vectors of doubles on the compiled machine, e.g.
    4 doubles if we use SSE2, 8 doubles for AVX. This
    is so that it matches up with MultiFloat, with both
    vectors providing the same number of elements.
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT MultiDouble
{

friend MultiDouble cos(const MultiDouble &val);
friend MultiDouble sin(const MultiDouble &val);
friend MultiDouble log(const MultiDouble &val);
friend MultiDouble exp(const MultiDouble &val);

friend void sincos(const MultiDouble &val, MultiDouble &sinval, MultiDouble &cosval);

public:
    MultiDouble();
    
    MultiDouble(double value);
    
    MultiDouble(const double *array, int size);
    MultiDouble(const QVector<float> &array);
    MultiDouble(const QVector<double> &array);
    
    MultiDouble(const MultiFloat &other);
    
    MultiDouble(const MultiDouble &other);
    
    #ifdef SIRE_HAS_CPP_11
        MultiDouble(const std::function<double ()> &func);
    #endif
    
    ~MultiDouble();
    
    bool isAligned() const;
    
    static QVector<MultiDouble> fromArray(const QVector<double> &array);
    static QVector<MultiDouble> fromArray(const QVector<float> &array);
    
    static QVector<MultiDouble> fromArray(const double *array, int size);
    static QVector<MultiDouble> fromArray(const float *array, int size);
    
    static QVector<double> toArray(const QVector<MultiDouble> &array);
    static QVector<double> toDoubleArray(const QVector<MultiDouble> &array);
    
    MultiDouble& operator=(const MultiDouble &other);
    MultiDouble& operator=(const MultiFloat &other);
    MultiDouble& operator=(double value);
    
    bool operator==(const MultiDouble &other) const;
    bool operator!=(const MultiDouble &other) const;
    
    bool operator<(const MultiDouble &other) const;
    bool operator>(const MultiDouble &other) const;
    
    bool operator<=(const MultiDouble &other) const;
    bool operator>=(const MultiDouble &other) const;
    
    MultiDouble compareEqual(const MultiDouble &other) const;
    MultiDouble compareNotEqual(const MultiDouble &other) const;

    MultiDouble compareLess(const MultiDouble &other) const;
    MultiDouble compareGreater(const MultiDouble &other) const;
    
    MultiDouble compareLessEqual(const MultiDouble &other) const;
    MultiDouble compareGreaterEqual(const MultiDouble &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    QString toBinaryString() const;
    
    static int size();
    static int count();
    
    double operator[](int i) const;
    
    void set(int i, double value);
    double get(int i) const;
    
    void quickSet(int i, double value);
    
    double at(int i) const;
    double getitem(int i) const;
    
    MultiDouble operator-() const;
    
    MultiDouble operator+(const MultiDouble &other) const;
    MultiDouble operator-(const MultiDouble &other) const;
    MultiDouble operator*(const MultiDouble &other) const;
    MultiDouble operator/(const MultiDouble &other) const;
    
    MultiDouble& operator+=(const MultiDouble &other);
    MultiDouble& operator-=(const MultiDouble &other);
    MultiDouble& operator*=(const MultiDouble &other);
    MultiDouble& operator/=(const MultiDouble &other);
    
    MultiDouble operator!() const;
    MultiDouble operator&(const MultiDouble &other) const;
    MultiDouble operator|(const MultiDouble &other) const;
    MultiDouble operator^(const MultiDouble &other) const;

    MultiDouble& operator&=(const MultiDouble &other);
    MultiDouble& operator|=(const MultiDouble &other);
    MultiDouble& operator^=(const MultiDouble &other);

    MultiDouble logicalNot() const;
    
    MultiDouble logicalAnd(const MultiDouble &other) const;
    MultiDouble logicalAndNot(const MultiDouble &other) const;
    
    MultiDouble logicalOr(const MultiDouble &other) const;
    MultiDouble logicalXor(const MultiDouble &other) const;
    
    MultiDouble& multiplyAdd(const MultiDouble &val0, const MultiDouble &val1);
    
    MultiDouble max(const MultiDouble &other) const;
    MultiDouble min(const MultiDouble &other) const;
    
    MultiDouble reciprocal() const;
    
    MultiDouble sqrt() const;
    MultiDouble rsqrt() const;
    MultiDouble rsqrt_approx() const;
    MultiDouble rsqrt_approx_nr() const;    

    MultiDouble rotate() const;
    
    double sum() const;
    double doubleSum() const;

    static void swap(MultiDouble &d0, int idx0, MultiDouble &d1, int idx1);

private:
    /* Give other Multi??? classes access to the raw vector */
    friend class MultiFloat;
    friend class MultiFixed;

    static void assertAligned(const void *ptr, size_t alignment);

    #ifndef SIRE_SKIP_INLINE_FUNCTIONS
        #ifndef MULTIFLOAT_CHECK_ALIGNMENT
            void assertAligned(){}
        #endif

        #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
            _ALIGNED(64) union
            {
                __m512d x[2];
                __m256d t[4];
                double a[16];
            } v;

            MultiDouble(__m512d avx_val0, __m512d avx_val1)
            {
                v.x[0] = avx_val0;
                v.x[1] = avx_val1;
            }

            MultiDouble(__mmask8 mask0, __mmask8 mask1)
            {
       	       	const __m512d zero = _mm512_set1_pd(0.0);
                const quint64 x = 0xFFFFFFFFFFFFFFFF;
       	       	const __m512d one = _mm512_set1_pd( 
                                *(reinterpret_cast<const double*>(&x)) );

                v.x[0] = _mm512_mask_blend_pd( mask0, zero, one );
                v.x[1] = _mm512_mask_blend_pd( mask1, zero, one );
            }

            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if ((quintptr)this % 64 != 0)
                        assertAligned(this, 64);
                }
            #endif
        #else
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            _ALIGNED(32) union
            {
                __m256d x[2];
                double a[8];
            } v;
        
            MultiDouble(__m256d avx_val0, __m256d avx_val1)
            {
                v.x[0] = avx_val0;
                v.x[1] = avx_val1;
            }
        
            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if ((quintptr)this % 32 != 0)
                        assertAligned(this, 32);
                }
            #endif
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            _ALIGNED(16) union
            {
                __m128d x[2];
                double a[4];
            } v;

            MultiDouble(__m128d sse_val0, __m128d sse_val1)
            {
                v.x[0] = sse_val0;
                v.x[1] = sse_val1;
            }

            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if ((quintptr)this % 16 != 0)
                        assertAligned(this, 16);
                }
            #endif
        #else
            _ALIGNED(32) union
            {
                double a[MULTIFLOAT_SIZE];
            } v;
            #define MULTIDOUBLE_BINONE getBinaryOne()

            static double getBinaryOne()
            {
                const quint64 x = 0xFFFFFFFFFFFFFFFFULL;
                return *(reinterpret_cast<const double*>(&x));
            }
        
            #ifdef MULTIFLOAT_CHECK_ALIGNMENT
                void assertAligned()
                {
                    if ((quintptr)this % 32 != 0)
                        assertAligned(this, 32);
                }
            #endif
        #endif
        #endif
        #endif
    #endif // #ifndef SIRE_SKIP_INLINE_FUNCTIONS
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

static const MultiDouble MULTIDOUBLE_ONE(1);

/** Constructor. This creates a MultiDouble with an undefined initial state */
inline
MultiDouble::MultiDouble()
{
    assertAligned();
}

/** Construct a vector with all values equal to 'val' */
inline
MultiDouble::MultiDouble(double val)
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_set1_pd(val);
        v.x[1] = _mm512_set1_pd(val);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_set1_pd(val);
        v.x[1] = _mm256_set1_pd(val);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_set1_pd(val);
        v.x[1] = _mm_set1_pd(val);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = val;
        }
    #endif
    #endif
    #endif
}

/** Construct using the passed function to generate doubles */
#ifdef SIRE_HAS_CPP_11
    inline
    MultiDouble::MultiDouble(const std::function<double ()> &generator)
    {
        assertAligned();

        #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
            v.x[0] = _mm512_set_pd(generator(), generator(), generator(), generator(),
                                   generator(), generator(), generator(), generator());
            v.x[1] = _mm512_set_pd(generator(), generator(), generator(), generator(),
                                   generator(), generator(), generator(), generator());
        #else
        #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
            v.x[0] = _mm256_set_pd(generator(), generator(), generator(), generator());
            v.x[1] = _mm256_set_pd(generator(), generator(), generator(), generator());
        #else
        #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
            v.x[0] = _mm_set_pd(generator(), generator());
            v.x[1] = _mm_set_pd(generator(), generator());
        #else
            for (int i=0; i<MULTIFLOAT_SIZE; ++i)
            {
                v.a[i] = generator();
            }
        #endif
        #endif
        #endif
    }
#endif

/** Copy construct from a MultiFloat */
inline
MultiDouble::MultiDouble(const MultiFloat &other)
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        const __m256 *o = (const __m256*)&(other.v.x);
        
        v.x[0] = _mm512_cvtps_pd( o[0] );
        v.x[1] = _mm512_cvtps_pd( o[1] );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        const __m128 *o = (const __m128*)&(other.v.x);
    
        v.x[0] = _mm256_cvtps_pd( o[0] );
        v.x[1] = _mm256_cvtps_pd( o[1] );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_cvtps_pd( other.v.x );
        v.x[1] = _mm_cvtps_pd( _mm_movehl_ps(other.v.x,other.v.x) );
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    #endif
    #endif
    #endif
}

/** Copy construct from a MultiDouble */
inline
MultiFloat::MultiFloat(const MultiDouble &other)
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        __m256 *o = (__m256*)&(v.x);

        o[0] = _mm512_cvtpd_ps(other.v.x[0]);
        o[1] = _mm512_cvtpd_ps(other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        __m128 *o = (__m128*)&(v.x);
    
        o[0] = _mm256_cvtpd_ps(other.v.x[0]);
        o[1] = _mm256_cvtpd_ps(other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_movelh_ps( _mm_cvtpd_ps(other.v.x[0]),
                             _mm_cvtpd_ps(other.v.x[1]) );
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    #endif
    #endif
    #endif
}

/** Copy constructor */
inline
MultiDouble::MultiDouble(const MultiDouble &other)
{
    assertAligned();

    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = other.v.x[0];
        v.x[1] = other.v.x[1];
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = other.v.x[0];
        v.x[1] = other.v.x[1];
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = other.v.x[0];
        v.x[1] = other.v.x[1];
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    #endif
    #endif
    #endif
}

/** Destructor */
inline
MultiDouble::~MultiDouble()
{}

/** Return the ith value in the MultiDouble - note that
    this is a fast function that does no bounds checking */
inline double MultiDouble::operator[](int i) const
{
    return v.a[i];
}

/** Assignment operator */
inline
MultiDouble& MultiDouble::operator=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = other.v.x[0];
        v.x[1] = other.v.x[1];
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = other.v.x[0];
        v.x[1] = other.v.x[1];
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = other.v.x[0];
        v.x[1] = other.v.x[1];
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    #endif
    #endif
    #endif
    
    return *this;
}

/** Assignment operator */
inline
MultiFloat& MultiFloat::operator=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        __m256 *o = (__m256*)&(v.x);

        o[0] = _mm512_cvtpd_ps(other.v.x[0]);
        o[1] = _mm512_cvtpd_ps(other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        __m128 *o = (__m128*)&(v.x);
    
        o[0] = _mm256_cvtpd_ps(other.v.x[0]);
        o[1] = _mm256_cvtpd_ps(other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x = _mm_movelh_ps( _mm_cvtpd_ps(other.v.x[0]),
                             _mm_cvtpd_ps(other.v.x[1]) );
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** Assignment operator */
inline
MultiDouble& MultiDouble::operator=(double value)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_set1_pd(value);
        v.x[1] = _mm512_set1_pd(value);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_set1_pd(value);
        v.x[1] = _mm256_set1_pd(value);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_set1_pd(value);
        v.x[1] = _mm_set1_pd(value);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = value;
        }
    #endif
    #endif
    #endif
    
    return *this;
}

/** Assignment operator */
inline
MultiDouble& MultiDouble::operator=(const MultiFloat &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        const __m256 *o = (const __m256*)&(other.v.x);

        v.x[0] = _mm512_cvtps_pd( o[0] );
        v.x[1] = _mm512_cvtps_pd( o[1] );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        const __m128 *o = (const __m128*)&(other.v.x);
    
        v.x[0] = _mm256_cvtps_pd( o[0] );
        v.x[1] = _mm256_cvtps_pd( o[1] );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_cvtps_pd( other.v.x );
        v.x[1] = _mm_cvtps_pd( _mm_movehl_ps(other.v.x,other.v.x) );
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] = other.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** Comparison operator. This will return a MultiDouble with elements
    set to zero for each double that is not equal */
inline
MultiDouble MultiDouble::compareEqual(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_cmp_pd_mask(v.x[0], other.v.x[0], _CMP_EQ_OQ),
                            _mm512_cmp_pd_mask(v.x[1], other.v.x[1], _CMP_EQ_OQ) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_cmp_pd(v.x[0], other.v.x[0], _CMP_EQ_OQ),
                            _mm256_cmp_pd(v.x[1], other.v.x[1], _CMP_EQ_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_cmpeq_pd(v.x[0], other.v.x[0]),
                            _mm_cmpeq_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;

        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] == other.v.a[i]) ? MULTIDOUBLE_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Not equals comparison operator */
inline
MultiDouble MultiDouble::compareNotEqual(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_cmp_pd_mask(v.x[0], other.v.x[0], _CMP_NEQ_OQ),
                            _mm512_cmp_pd_mask(v.x[1], other.v.x[1], _CMP_NEQ_OQ) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_cmp_pd(v.x[0], other.v.x[0], _CMP_NEQ_OQ),
                            _mm256_cmp_pd(v.x[1], other.v.x[1], _CMP_NEQ_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_cmpneq_pd(v.x[0], other.v.x[0]),
                            _mm_cmpneq_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] != other.v.a[i]) ? MULTIDOUBLE_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Less than comparison operator */
inline
MultiDouble MultiDouble::compareLess(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_cmp_pd_mask(v.x[0], other.v.x[0], _CMP_LT_OQ),
                            _mm512_cmp_pd_mask(v.x[1], other.v.x[1], _CMP_LT_OQ) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_cmp_pd(v.x[0], other.v.x[0], _CMP_LT_OQ),
                            _mm256_cmp_pd(v.x[1], other.v.x[1], _CMP_LT_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_cmplt_pd(v.x[0], other.v.x[0]),
                            _mm_cmplt_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] < other.v.a[i]) ? MULTIDOUBLE_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Greater than comparison operator */
inline
MultiDouble MultiDouble::compareGreater(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_cmp_pd_mask(v.x[0], other.v.x[0], _CMP_GT_OQ),
                            _mm512_cmp_pd_mask(v.x[1], other.v.x[1], _CMP_GT_OQ) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_cmp_pd(v.x[0], other.v.x[0], _CMP_GT_OQ),
                            _mm256_cmp_pd(v.x[1], other.v.x[1], _CMP_GT_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_cmpgt_pd(v.x[0], other.v.x[0]),
                            _mm_cmpgt_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] > other.v.a[i]) ? MULTIDOUBLE_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Less than or equal comparison */
inline
MultiDouble MultiDouble::compareLessEqual(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_cmp_pd_mask(v.x[0], other.v.x[0], _CMP_LE_OQ),
                            _mm512_cmp_pd_mask(v.x[1], other.v.x[1], _CMP_LE_OQ) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_cmp_pd(v.x[0], other.v.x[0], _CMP_LE_OQ),
                            _mm256_cmp_pd(v.x[1], other.v.x[1], _CMP_LE_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_cmple_pd(v.x[0], other.v.x[0]),
                            _mm_cmple_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] <= other.v.a[i]) ? MULTIDOUBLE_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Greater than or equal comparison */
inline
MultiDouble MultiDouble::compareGreaterEqual(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_cmp_pd_mask(v.x[0], other.v.x[0], _CMP_GE_OQ),
                            _mm512_cmp_pd_mask(v.x[1], other.v.x[1], _CMP_GE_OQ) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_cmp_pd(v.x[0], other.v.x[0], _CMP_GE_OQ),
                            _mm256_cmp_pd(v.x[1], other.v.x[1], _CMP_GE_OQ) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_cmpge_pd(v.x[0], other.v.x[0]),
                            _mm_cmpge_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = (v.a[i] >= other.v.a[i]) ? MULTIDOUBLE_BINONE : 0x0;
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Return the number of values in the vector */
inline
int MultiDouble::size()
{
    return MULTIFLOAT_SIZE;
}

/** Return the number of values in the vector */
inline
int MultiDouble::count()
{
    return MULTIFLOAT_SIZE;
}

/** Addition operator */
inline
MultiDouble MultiDouble::operator+(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_add_pd(v.x[0], other.v.x[0]),
                            _mm512_add_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_add_pd(v.x[0], other.v.x[0]),
                            _mm256_add_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_add_pd(v.x[0], other.v.x[0]),
                            _mm_add_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] + other.v.a[i];
        }
        return ret;
    #endif
    #endif
    #endif
}

/** Subtraction operator */
inline
MultiDouble MultiDouble::operator-(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_sub_pd(v.x[0], other.v.x[0]),
                            _mm512_sub_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_sub_pd(v.x[0], other.v.x[0]),
                            _mm256_sub_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_sub_pd(v.x[0], other.v.x[0]),
                            _mm_sub_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] - other.v.a[i];
        }
        return ret;
    #endif
    #endif
    #endif
}

/** Multiplication operator */
inline
MultiDouble MultiDouble::operator*(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_mul_pd(v.x[0], other.v.x[0]),
                            _mm512_mul_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_mul_pd(v.x[0], other.v.x[0]),
                            _mm256_mul_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_mul_pd(v.x[0], other.v.x[0]),
                            _mm_mul_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] * other.v.a[i];
        }
        return ret;
    #endif
    #endif
    #endif
}

/** Division operator */
inline
MultiDouble MultiDouble::operator/(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_div_pd(v.x[0], other.v.x[0]),
                            _mm512_div_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_div_pd(v.x[0], other.v.x[0]),
                            _mm256_div_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_div_pd(v.x[0], other.v.x[0]),
                            _mm_div_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = v.a[i] / other.v.a[i];
        }
        return ret;
    #endif
    #endif
    #endif
}

/** In-place addition operator */
inline
MultiDouble& MultiDouble::operator+=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_add_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm512_add_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_add_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_add_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_add_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_add_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] += other.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** In-place subtraction operator */
inline
MultiDouble& MultiDouble::operator-=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_sub_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm512_sub_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_sub_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_sub_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_sub_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_sub_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] -= other.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** In-place multiplication operator */
inline
MultiDouble& MultiDouble::operator*=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_mul_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm512_mul_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_mul_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_mul_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_mul_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_mul_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] *= other.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** In-place division operator */
inline
MultiDouble& MultiDouble::operator/=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_div_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm512_div_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_div_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_div_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_div_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_div_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] /= other.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** Bitwise logical "and" comparison */
inline
MultiDouble MultiDouble::logicalAnd(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        MultiDouble ret;
        ret.v.t[0] = _mm256_and_pd(v.t[0], other.v.t[0]);
        ret.v.t[1] = _mm256_and_pd(v.t[1], other.v.t[1]);
        ret.v.t[2] = _mm256_and_pd(v.t[2], other.v.t[2]);
        ret.v.t[3] = _mm256_and_pd(v.t[3], other.v.t[3]);
        return ret;
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_and_pd(v.x[0], other.v.x[0]),
                            _mm256_and_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_and_pd(v.x[0], other.v.x[0]),
                            _mm_and_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                ret_char_v[j] = char_v[j] & other_char_v[j];
            }
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Bitwise logical "and not" */
inline
MultiDouble MultiDouble::logicalAndNot(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        MultiDouble ret;
        ret.v.t[0] = _mm256_andnot_pd(v.t[0], other.v.t[0]);
        ret.v.t[1] = _mm256_andnot_pd(v.t[1], other.v.t[1]);
        ret.v.t[2] = _mm256_andnot_pd(v.t[2], other.v.t[2]);
        ret.v.t[3] = _mm256_andnot_pd(v.t[3], other.v.t[3]);
        return ret;
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_andnot_pd(v.x[0], other.v.x[0]),
                            _mm256_andnot_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_andnot_pd(v.x[0], other.v.x[0]),
                            _mm_andnot_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                ret_char_v[j] = !(char_v[j] & other_char_v[j]);
            }
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Bitwise logical or operator */
inline
MultiDouble MultiDouble::logicalOr(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_or_pd(v.x[0], other.v.x[0]),
                            _mm512_or_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_or_pd(v.x[0], other.v.x[0]),
                            _mm256_or_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_or_pd(v.x[0], other.v.x[0]),
                            _mm_or_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                ret_char_v[j] = char_v[j] | other_char_v[j];
            }
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Bitwise logical xor */
inline
MultiDouble MultiDouble::logicalXor(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_xor_pd(v.x[0], other.v.x[0]),
                            _mm512_xor_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_xor_pd(v.x[0], other.v.x[0]),
                            _mm256_xor_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_xor_pd(v.x[0], other.v.x[0]),
                            _mm_xor_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
    
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
            const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                ret_char_v[j] = char_v[j] ^ other_char_v[j];
            }
        }
    
        return ret;
    #endif
    #endif
    #endif
}

/** Logical not operator */
inline
MultiDouble MultiDouble::logicalNot() const
{
    MultiDouble ret;

    for (int i=0; i<MULTIFLOAT_SIZE; ++i)
    {
        unsigned char *ret_char_v = reinterpret_cast<unsigned char*>(&(ret.v.a[i]));
        const unsigned char *char_v = reinterpret_cast<const unsigned char*>(&(v.a[i]));

        for (unsigned int j=0; j<sizeof(double); ++j)
        {
            ret_char_v[j] = !char_v[j];
        }
    }

    return ret;
}

/** Logical not operator */
inline
MultiDouble MultiDouble::operator!() const
{
    return this->logicalNot();
}

/** Logical and operator */
inline
MultiDouble MultiDouble::operator&(const MultiDouble &other) const
{
    return this->logicalAnd(other);
}

/** Logical or operator */
inline
MultiDouble MultiDouble::operator|(const MultiDouble &other) const
{
    return this->logicalOr(other);
}

/** Logical xor operator */
inline
MultiDouble MultiDouble::operator^(const MultiDouble &other) const
{
    return this->logicalXor(other);
}

/** In place logical and */
inline
MultiDouble& MultiDouble::operator&=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.t[0] = _mm256_and_pd(v.t[0], other.v.t[0]);
        v.t[1] = _mm256_and_pd(v.t[1], other.v.t[1]);
       	v.t[2] = _mm256_and_pd(v.t[2], other.v.t[2]);
       	v.t[3] = _mm256_and_pd(v.t[3], other.v.t[3]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_and_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_and_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_and_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_and_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                char_v[j] &= other_char_v[j];
            }
        }
    #endif
    #endif
    #endif

    return *this;
}

/** In-place logical or */
inline
MultiDouble& MultiDouble::operator|=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_or_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm512_or_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_or_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_or_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_or_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_or_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                char_v[j] |= other_char_v[j];
            }
        }
    #endif
    #endif
    #endif

    return *this;
}

/** In-place logical xor */
inline
MultiDouble& MultiDouble::operator^=(const MultiDouble &other)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_xor_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm512_xor_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_xor_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm256_xor_pd(v.x[1], other.v.x[1]);
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_xor_pd(v.x[0], other.v.x[0]);
        v.x[1] = _mm_xor_pd(v.x[1], other.v.x[1]);
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            unsigned char *char_v = reinterpret_cast<unsigned char*>(&(v.a[i]));
            const unsigned char *other_char_v
                        = reinterpret_cast<const unsigned char*>(&(other.v.a[i]));

            for (unsigned int j=0; j<sizeof(double); ++j)
            {
                char_v[j] ^= other_char_v[j];
            }
        }
    #endif
    #endif
    #endif

    return *this;
}

/** Multiply val0 and val1 and add it onto this value */
inline
MultiDouble& MultiDouble::multiplyAdd(const MultiDouble &v0, const MultiDouble &v1)
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        v.x[0] = _mm512_add_pd(v.x[0], _mm512_mul_pd(v0.v.x[0], v1.v.x[0]));
        v.x[1] = _mm512_add_pd(v.x[1], _mm512_mul_pd(v0.v.x[1], v1.v.x[1]));
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        v.x[0] = _mm256_add_pd(v.x[0], _mm256_mul_pd(v0.v.x[0], v1.v.x[0]));
        v.x[1] = _mm256_add_pd(v.x[1], _mm256_mul_pd(v0.v.x[1], v1.v.x[1]));
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        v.x[0] = _mm_add_pd(v.x[0], _mm_mul_pd(v0.v.x[0], v1.v.x[0]));
        v.x[1] = _mm_add_pd(v.x[1], _mm_mul_pd(v0.v.x[1], v1.v.x[1]));
    #else
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            v.a[i] += v0.v.a[i] * v1.v.a[i];
        }
    #endif
    #endif
    #endif

    return *this;
}

/** Return the maximum vector between this and other */
inline
MultiDouble MultiDouble::max(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_max_pd(v.x[0], other.v.x[0]),
                            _mm512_max_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_max_pd(v.x[0], other.v.x[0]),
                            _mm256_max_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_max_pd(v.x[0], other.v.x[0]),
                            _mm_max_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::max(v.a[i], other.v.a[i]);
        }
        return ret;
    #endif
    #endif
    #endif
}

/** Return the minimum vector between this and other */
inline
MultiDouble MultiDouble::min(const MultiDouble &other) const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_min_pd(v.x[0], other.v.x[0]),
                            _mm512_min_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_min_pd(v.x[0], other.v.x[0]),
                            _mm256_min_pd(v.x[1], other.v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_min_pd(v.x[0], other.v.x[0]),
                            _mm_min_pd(v.x[1], other.v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::min(v.a[i], other.v.a[i]);
        }
        return ret;
    #endif
    #endif
    #endif
}

/** Return the reciprocal of this vector */
inline
MultiDouble MultiDouble::reciprocal() const
{
    return MULTIDOUBLE_ONE.operator/(*this);
}

/** Return the square root of this vector */
inline
MultiDouble MultiDouble::sqrt() const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_sqrt_pd(v.x[0]),
                            _mm512_sqrt_pd(v.x[1]) );
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return MultiDouble( _mm256_sqrt_pd(v.x[0]),
                            _mm256_sqrt_pd(v.x[1]) );
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return MultiDouble( _mm_sqrt_pd(v.x[0]),
                            _mm_sqrt_pd(v.x[1]) );
    #else
        MultiDouble ret;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i] = std::sqrt(v.a[i]);
        }
        return ret;
    #endif
    #endif
    #endif
}

/** Return the recipical square root of this vector */
inline
MultiDouble MultiDouble::rsqrt() const
{
    return MULTIDOUBLE_ONE.operator/(this->sqrt());
}

/** Return an approximation of the reciprical square root of this vector (only good
    for about 12 bit of the mantissa) */
inline
MultiDouble MultiDouble::rsqrt_approx() const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return MultiDouble( _mm512_rsqrt14_pd(v.x[0]), _mm512_rsqrt14_pd(v.x[1]) );
    #else
        return MULTIFLOAT_ONE.operator/(this->sqrt());
    #endif
}

/** Return a good approximation of the reciprical square root (this poor approximation
    refined using a single Newton Raphson step) */
inline
MultiDouble MultiDouble::rsqrt_approx_nr() const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        //get the approximation
        __m512d a0 = _mm512_rsqrt14_pd(v.x[0]);

        //now use one step of NR to refine the result
        // 1/x = 0.5 a[ 3 - x a^2 ] where a is the approximation

        __m512d tmp = _mm512_mul_pd(a0, v.x[0]);
        tmp = _mm512_mul_pd(a0, tmp);
        const __m512d three = _mm512_set1_pd(3.0);
        tmp = _mm512_sub_pd(three, tmp);
        a0 = _mm512_mul_pd(a0, tmp);
        const __m512d half = _mm512_set1_pd(0.5);
        a0 = _mm512_mul_pd(a0, half);

        //repeat for the other double
        __m512d a1 = _mm512_rsqrt14_pd(v.x[1]);
        tmp = _mm512_mul_pd(a1, v.x[1]);
        tmp = _mm512_mul_pd(a1, tmp); 
        tmp = _mm512_sub_pd(three, tmp);
        a1 = _mm512_mul_pd(a1, tmp);
        a1 = _mm512_mul_pd(a1, half);

        return MultiDouble(a0, a1);
    #else
        return this->rsqrt();
    #endif
}

/** Rotate this vector. This moves each element one space to the left, moving the
    first element to the last element */
inline
MultiDouble MultiDouble::rotate() const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        MultiDouble ret;

        ret.v.x[0] = _mm512_permutexvar_pd(
           _mm512_set_epi32(0, 0, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1),
               v.x[0] );

        ret.v.x[1] = _mm512_permutexvar_pd(
           _mm512_set_epi32(0, 0, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1),
               v.x[1] );

        ret.v.a[15] = v.a[0];
        ret.v.a[7] = v.a[8];

        return ret;
    #else
        MultiDouble ret;
    
        for (int i=1; i<MULTIFLOAT_SIZE; ++i)
        {
            ret.v.a[i-1] = v.a[i];
        }
    
        ret.v.a[MULTIFLOAT_SIZE-1] = v.a[0];
 
        return ret;
    #endif
}

/** Return the sum of all elements of this vector */
inline
double MultiDouble::sum() const
{
    #ifdef MULTIFLOAT_AVX512F_IS_AVAILABLE
        return v.a[0] + v.a[1] + v.a[2] + v.a[3] +
               v.a[4] + v.a[5] + v.a[6] + v.a[7] +
               v.a[8] + v.a[9] + v.a[10] + v.a[11] +
               v.a[12] + v.a[13] + v.a[14] + v.a[15];
    #else
    #ifdef MULTIFLOAT_AVX_IS_AVAILABLE
        return v.a[0] + v.a[1] + v.a[2] + v.a[3] +
               v.a[4] + v.a[5] + v.a[6] + v.a[7];
    #else
    #ifdef MULTIFLOAT_SSE_IS_AVAILABLE
        return v.a[0] + v.a[1] + v.a[2] + v.a[3];
    #else
        double sum = 0;
        for (int i=0; i<MULTIFLOAT_SIZE; ++i)
        {
            sum += v.a[i];
        }
        return sum;
    #endif
    #endif
    #endif
}

/** Return the sum of all elements of this vector, using doubles for the sum */
inline
double MultiDouble::doubleSum() const
{
    return this->sum();
}

inline MultiDouble cos(const MultiDouble &val)
{
    return cos(MultiFloat(val));
}

inline MultiDouble sin(const MultiDouble &val)
{
    return sin(MultiFloat(val));
}

inline MultiDouble log(const MultiDouble &val)
{
    return log(MultiFloat(val));
}

inline MultiDouble exp(const MultiDouble &val)
{
    return exp(MultiFloat(val));
}

inline void sincos(const MultiDouble &val, MultiDouble &sinval, MultiDouble &cosval)
{
    MultiFloat c, s;
    sincos(MultiFloat(val), c, s);
    sinval = s;
    cosval = c;
}

#endif // #ifndef SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireMaths::MultiDouble )

SIRE_END_HEADER

#endif

