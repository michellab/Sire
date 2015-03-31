/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMATHS_MATHS_H
#define SIREMATHS_MATHS_H

#include "sireglobal.h"

#include <cmath>

#ifndef GCCXML_PARSE
    //some versions of GCCXML don't like these headers...
    #include <boost/random/variate_generator.hpp>
    #include <boost/random/uniform_real.hpp>
    #include <boost/random/uniform_int.hpp>
    #include <boost/random/mersenne_twister.hpp>
#endif // GCCXML_PARSE

#include <gsl/gsl_sys.h>

#include "constants.h"

#ifdef SIRE_MAX
#undef SIRE_MAX
#endif

#ifdef SIRE_MIN
#undef SIRE_MIN
#endif

#define SIRE_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define SIRE_MIN(x, y) (((x) < (y)) ? (x) : (y))

SIRE_BEGIN_HEADER

namespace SireMaths
{

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

typedef boost::variate_generator<boost::mt19937, boost::uniform_real<double> > UniformRand;
typedef boost::variate_generator<boost::mt19937, boost::uniform_int<qint32> > UniformRanInt;
typedef boost::variate_generator<boost::mt19937, boost::uniform_int<quint32> > UniformRanUInt;
typedef boost::variate_generator<boost::mt19937, boost::uniform_int<qint64> > UniformRanInt64;
typedef boost::variate_generator<boost::mt19937, boost::uniform_int<quint64> > UniformRanUInt64;

/** Function to return a uniform double random number generator with a range
    from 'min' to 'max' */
inline UniformRand uniformRanGenerator(double min=0.0, double max=1.0)
{
    return UniformRand( boost::mt19937(), boost::uniform_real<double>(min,max) );
}

/** Function to return a uniform double random number generator with a range
    from 'min' to 'max' and seeded with 'seed' */
inline UniformRand uniformRanGenerator(double min, double max, quint32 seed)
{
    return UniformRand( boost::mt19937(seed), boost::uniform_real<double>(min,max) );
}

/** Function to return a uniform int random number generate with a range from
    'min' to 'max' */
inline UniformRanInt uniformRanIntGenerator(qint32 min=0, qint32 max=100)
{
    return UniformRanInt( boost::mt19937(), boost::uniform_int<qint32>(min,max) );
}

/** Function to return a uniform int random number generate with a range from
    'min' to 'max' and seeded with 'seed' */
inline UniformRanInt uniformRanIntGenerator(qint32 min, qint32 max, quint32 seed)
{
    return UniformRanInt( boost::mt19937(seed), boost::uniform_int<qint32>(min,max) );
}

/** Function to return a uniform uint random number generate with a range from
    'min' to 'max' */
inline UniformRanUInt uniformRanUIntGenerator(quint32 min=0, quint32 max=100)
{
    return UniformRanUInt( boost::mt19937(), boost::uniform_int<quint32>(min,max) );
}

/** Function to return a uniform uint random number generate with a range from
    'min' to 'max' and seeded with 'seed' */
inline UniformRanUInt uniformRanUIntGenerator(quint32 min, quint32 max, quint32 seed)
{
    return UniformRanUInt( boost::mt19937(seed), boost::uniform_int<quint32>(min,max) );
}

/** Function to return a uniform int random number generate with a range from
    'min' to 'max' */
inline UniformRanInt64 uniformRanInt64Generator(qint64 min=0, qint64 max=100)
{
    return UniformRanInt64( boost::mt19937(), boost::uniform_int<qint64>(min,max) );
}

/** Function to return a uniform int random number generate with a range from
    'min' to 'max' and seeded with 'seed' */
inline UniformRanInt64 uniformRanInt64Generator(qint64 min, qint64 max, quint32 seed)
{
    return UniformRanInt64( boost::mt19937(seed), boost::uniform_int<qint64>(min,max) );
}

/** Function to return a uniform uint random number generate with a range from
    'min' to 'max' */
inline UniformRanUInt64 uniformRanUInt64Generator(quint64 min=0, quint64 max=100)
{
    return UniformRanUInt64( boost::mt19937(), boost::uniform_int<quint64>(min,max) );
}

/** Function to return a uniform uint random number generate with a range from
    'min' to 'max' and seeded with 'seed' */
inline UniformRanUInt64 uniformRanUInt64Generator(quint64 min, quint64 max, quint32 seed)
{
    return UniformRanUInt64( boost::mt19937(seed), boost::uniform_int<quint64>(min,max) );
}

/** Return x raised to the power 2 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_2(const T &x) { return x*x;   }

/** Return x raised to the power 3 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_3(const T &x) { return x*x*x; }

/** Return x raised to the power 4 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_4(const T &x) { T x2 = x*x;   return x2*x2;    }

/** Return x raised to the power 5 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_5(const T &x) { T x2 = x*x;   return x2*x2*x;  }

/** Return x raised to the power 6 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_6(const T &x) { T x2 = x*x;   return x2*x2*x2; }

/** Return x raised to the power 7 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_7(const T &x) { T x3 = x*x*x; return x3*x3*x;  }

/** Return x raised to the power 8 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_8(const T &x) { T x2 = x*x;   T x4 = x2*x2; return x4*x4; }

/** Return x raised to the power 9 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_9(const T &x) { T x3 = x*x*x; return x3*x3*x3; }

/** Return x raised to the power 10 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_10(const T &x) { T x2 = x*x; T x5 = x2*x2*x; return x5*x5; }

/** Return x raised to the power 11 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_11(const T &x) { T x2 = x*x; T x5 = x2*x2*x; return x5*x5*x; }

/** Return x raised to the power 12 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_12(const T &x) { T x2 = x*x; T x6 = x2*x2*x2; return x6*x6; }

/** Return x raised to the power 13 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_13(const T &x) { T x2 = x*x; T x6 = x2*x2*x2; return x6*x6*x; }

/** Return x raised to the power 14 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_14(const T &x) { T x2 = x*x; T x7 = x2*x2*x2*x; return x7*x7; }

/** Return x raised to the power 15 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_15(const T &x) { T x2 = x*x; T x5 = x2*x2*x; return x5*x5*x5; }

/** Return x raised to the power 16 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_16(const T &x) { T x2 = x*x; T x4 = x2*x2; return x4*x4*x4*x4; }

/** Return x raised to the power 17 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_17(const T &x) { T x2 = x*x; T x4 = x2*x2; return x4*x4*x4*x4*x; }

/** Return x raised to the power 18 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_18(const T &x) { T x3 = x*x*x; T x9 = x3*x3*x3; return x9*x9; }

/** Return x raised to the power 19 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_19(const T &x) { T x3 = x*x*x; T x9 = x3*x3*x3; return x9*x9*x; }

/** Return x raised to the power 20 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_20(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x8=x4*x4; return x8*x8*x4; }

/** Return x raised to the power 21 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_21(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x7=x4*x2*x; return x7*x7*x7; }

/** Return x raised to the power 22 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_22(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x6=x4*x2; return x6*x6*x4; }

/** Return x raised to the power 23 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_23(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x6=x4*x2; return x6*x6*x4*x; }

/** Return x raised to the power 24 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_24(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x8=x4*x4; return x8*x8*x8; }

/** Internal function used by the 'pow' function - raise 'x' to
    the positive power 'n' */
inline double pow_pvt(double x, int n)
{
    switch(n)
    {
        case 1:
            return x;
        case 2:
            return pow_2(x);
        case 3:
            return pow_3(x);
        case 4:
            return pow_4(x);
        case 5:
            return pow_5(x);
        case 6:
            return pow_6(x);
        case 7:
            return pow_7(x);
        case 8:
            return pow_8(x);
        case 9:
            return pow_9(x);
        case 10:
            return pow_10(x);
        case 11:
            return pow_11(x);
        case 12:
            return pow_12(x);
        default:
            return std::pow(x, n);
    }
}

/** Return x raised to the power n */
inline double pow(double x, int n)
{
    if (n > 0)
        return pow_pvt(x,n);
    else if (n == 0)
        return double(1.0);
    else
        return double(1.0) / pow_pvt(x, -1*n);
}

/** Return whether or not an integer is even */
inline bool isEven(int val)
{
    return val % 2 == 0;
}

/** Return whether or not an integer is odd */
inline bool isOdd(int val)
{
    return val % 2 != 0;
}

/** Return whether or not this is an integer */
inline bool isInteger(double val)
{
    return (smallest_int <= val) and (val <= largest_int) and (val == double(int(val)));
}

/** Return pow(x,n) for non-integer n */
inline double pow(double x, double n)
{
    if (isInteger(n))
        return pow(x, int(n));
    else
        return std::pow(x,n);
}

/** Here is a non-fast version of fastInvSqrt */
inline double invSqrt(double x)
{
    return double(1.0) / std::sqrt(x);
}

/** Here is a non-fast version of fastInvSqrt */
inline float invSqrt(float x)
{
    return float(1.0) / std::sqrt(x);
}

/** Return true if two numbers are equal. This uses gsl_fcmp
    for the comparison, and the values must be equal to within
    a range of 2*delta, where delta = 2^k * epsilon, where
    k equals the maximum base-2 exponent of val0 and val1 as
    calculated via frexp(). This should normally be a value
    around 1e-6 */
inline bool areEqual(double val0, double val1, double epsilon)
{
    return not ( gsl_fcmp(val0,val1,epsilon) );
}

/** Return true if two numbers are equal */
inline bool areEqual(double val0, double val1)
{
    return not ( gsl_fcmp(val0,val1,1e-6) );
}

/** Return true if this number is equal to zero */
template<typename T>
SIRE_INLINE_TEMPLATE
bool isZero(const T &val)
{
    return (std::abs(val) < std::numeric_limits<T>::epsilon());
}

#else  // SIRE_SKIP_INLINE_FUNCTIONS

class UniformRand;
class UniformRanInt;
class UniformRanUInt;
class UniformRanInt64;
class UniformRanUInt64;

UniformRand uniformRanGenerator(double min=0.0, double max=1.0);
UniformRand uniformRanGenerator(double min, double max, quint32 seed);
UniformRanInt uniformRanIntGenerator(qint32 min=0, qint32 max=100);
UniformRanInt uniformRanIntGenerator(qint32 min, qint32 max, quint32 seed);
UniformRanUInt uniformRanUIntGenerator(quint32 min=0, quint32 max=100);
UniformRanUInt uniformRanUIntGenerator(quint32 min, quint32 max, quint32 seed);
UniformRanInt64 uniformRanInt64Generator(qint64 min=0, qint64 max=100);
UniformRanInt64 uniformRanInt64Generator(qint64 min, qint64 max, quint32 seed);
UniformRanUInt64 uniformRanUInt64Generator(quint64 min=0, quint64 max=100);
UniformRanUInt64 uniformRanUInt64Generator(quint64 min, quint64 max, quint32 seed);

template<typename T>
T pow_2(const T &x);

template<typename T>
T pow_3(const T &x);

template<typename T>
T pow_4(const T &x);

template<typename T>
T pow_5(const T &x);

template<typename T>
T pow_6(const T &x);

template<typename T>
T pow_7(const T &x);

template<typename T>
T pow_8(const T &x);

template<typename T>
T pow_9(const T &x);

template<typename T>
T pow_10(const T &x);

template<typename T>
T pow_11(const T &x);

template<typename T>
T pow_12(const T &x);

template<typename T>
T pow_13(const T &x);

template<typename T>
T pow_14(const T &x);

template<typename T>
T pow_15(const T &x);

template<typename T>
T pow_16(const T &x);

template<typename T>
T pow_17(const T &x);

template<typename T>
T pow_18(const T &x);

template<typename T>
T pow_19(const T &x);

template<typename T>
T pow_20(const T &x);

template<typename T>
T pow_21(const T &x);

template<typename T>
T pow_22(const T &x);

template<typename T>
T pow_23(const T &x);

template<typename T>
T pow_24(const T &x);

double pow_pvt(double x, int n);

double pow(double x, int n);

bool isEven(int val);
bool isOdd(int val);
bool isInteger(double val);
double pow(double x, double n);

double invSqrt(double x);
float invSqrt(float x);

bool areEqual(double val0, double val1, double epsilon);
bool areEqual(double val0, double val1);

template<typename T>
bool isZero(const T &val);

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
