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

#ifndef SIREMATHS_VECTOR_H
#define SIREMATHS_VECTOR_H

class QDataStream;
class QString;

#include <QObject>

#include <boost/tuple/tuple.hpp>

#include "SireMaths/errors.h"

#include "maths.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Vector;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::Vector&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::Vector&);

uint qHash(const SireMaths::Vector&);

namespace SireMaths
{

class Quaternion;
class Matrix;
class NVector;

using boost::tuple;

using SireUnits::Dimension::Angle;
using SireUnits::Dimension::Length;

SIREMATHS_EXPORT const Vector operator+(const Vector &p1, const Vector &p2);
SIREMATHS_EXPORT const Vector operator-(const Vector &p1, const Vector &p2);
SIREMATHS_EXPORT const Vector operator*(const Vector &p1, double c);
SIREMATHS_EXPORT const Vector operator*(double c, const Vector &p1);
SIREMATHS_EXPORT const Vector operator/(const Vector &p1, double c);
SIREMATHS_EXPORT const Quaternion operator*(const Vector &p1, const Quaternion &p2);
SIREMATHS_EXPORT const Quaternion operator*(const Quaternion &p1, const Vector &p2);
SIREMATHS_EXPORT const Vector operator*(const Matrix &m, const Vector &p);

/**
This is a simple 3D vector.

A vector is implemented as a very simple collection of 3 doubles. It should
thus be very efficient in terms of speed and size as there is no array. This
does mean that the .x(), .y() and .z() functions provide faster access
than the [] operator.

There is a length/distance and a fastLength/fastDistance. The normal
version uses the standard library's sqrt function, while the fast versions
use a numerical approximation that is quite accurate (to within 1-2%).

These have been benchmarked on my development system using the 'testSpeed'
static function. The results are;

100,000,000 evaluations of vectors (1.0,2.0,3.0), (3.0,2.0,1.0)

distance took on average 4498 ms
distance2 took on average 148 ms
fastDistance took on average 350 ms

This shows that you should normally use 'distance2' for distance checking,
with fastDistance if you are not concerned with total accuracy.

The error with fastDistance has an oscillating value, with the error having
a maximum overestimate of 1*10-5, and maximum underestimate of 0.00065 between
the distances of 0.0 and 150.0

Within this range the percentage error also oscillates, with the oscillations
having a random frequency, but almost constant amplitude of 0.0005%. See the
graph in techdocs/fastLengthError.png and the script techdocs/fastLengthError.py

@author Christopher Woods
*/

class Quaternion;

class SIREMATHS_EXPORT Vector
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const Vector&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, Vector&);

public:
    typedef double value_type;

    Vector( double val=0.0 );
    Vector( Length val );
    Vector( double xpos, double ypos, double zpos );
    Vector( Length xpos, Length ypos, Length zpos );
    Vector( const tuple<double,double,double> &pos );
    Vector( const tuple<Length,Length,Length> &pos );

    Vector( const QString &str );

    Vector(const NVector &other);

    Vector(const Vector &other);

    ~Vector();

    static const char* typeName();

    const char* what() const
    {
        return Vector::typeName();
    }

    double* data();

    const double* data() const;
    const double* constData() const;

    double x() const;
    double y() const;
    double z() const;

    double r() const;
    double g() const;
    double b() const;

    Vector direction() const;
    double magnitude() const;

    const Vector& operator=(const Vector &other);

    bool operator==(const Vector &p1) const;
    bool operator!=(const Vector &p1) const;

    const Vector& operator+=(const Vector &other);
    const Vector& operator-=(const Vector &other);
    const Vector& operator*=(const double &other);
    const Vector& operator/=(const double &other);
    Vector operator-() const;

    void set(double x, double y, double z);
    void setX(double x);
    void setY(double y);
    void setZ(double z);
    void setR(double x);
    void setG(double y);
    void setB(double z);

    void set(Length x, Length y, Length z);
    void setX(Length x);
    void setY(Length y);
    void setZ(Length z);

    void set(int i, double v);
    void set(int i, Length v);

    double operator[](int i) const;

    int count() const;

    double at(int i) const;

    double getitem(int i) const;

    double manhattanLength() const;

    double length() const;
    double length2() const;

    double invLength() const;
    double invLength2() const;

    Vector normalise() const;

    bool isZero() const;

    QString toString() const;

    static Vector fromString(const QString &str);

    static double dot(const Vector &v0, const Vector &v1);
    // Note that the cross product returns a normal vector, i.e. it has been
    // normalised!
    static Vector cross(const Vector &v0, const Vector &v1);

    // This is a regular cross product.
    static Vector realCross(const Vector &v0, const Vector &v1);

    void setMax(const Vector &other);
    void setMin(const Vector &other);

    Vector max(const Vector &other) const;
    Vector min(const Vector &other) const;

    Angle bearing() const;
    Angle bearingXY(const Vector &v) const;
    Angle bearingXZ(const Vector &v) const;
    Angle bearingYZ(const Vector &v) const;

    Matrix metricTensor() const;

    static double distance2(const Vector &v1, const Vector &v2);
    static double distance(const Vector &v1, const Vector &v2);

    static double invDistance(const Vector &v1, const Vector &v2);
    static double invDistance2(const Vector &v1, const Vector &v2);

    static Angle angle(const Vector &v0, const Vector &v1);
    static Angle angle(const Vector &v0, const Vector &v1, const Vector &v2);

    static Angle dihedral(const Vector &v0, const Vector &v1,
                          const Vector &v2, const Vector &v3);

    static Vector generate(double dst, const Vector &v1, const Angle &ang,
                           const Vector &v2, const Angle &dih, const Vector &v3);

    static Vector generate(Length dst, const Vector &v1, const Angle &ang,
                           const Vector &v2, const Angle &dih, const Vector &v3);

    friend SIREMATHS_EXPORT const Vector operator+(const Vector &p1, const Vector &p2);
    friend SIREMATHS_EXPORT const Vector operator-(const Vector &p1, const Vector &p2);
    friend SIREMATHS_EXPORT const Vector operator*(const Vector &p1, double c);
    friend SIREMATHS_EXPORT const Vector operator*(double c, const Vector &p1);
    friend SIREMATHS_EXPORT const Vector operator/(const Vector &p1, double c);
    friend SIREMATHS_EXPORT const Quaternion operator*(const Vector &p1, const Quaternion &p2);
    friend SIREMATHS_EXPORT const Quaternion operator*(const Quaternion &p1, const Vector &p2);
    friend SIREMATHS_EXPORT const Vector operator*(const Matrix &m, const Vector &p);

protected:
    /** Use four values so that vector arrays can be nicely aligned */
    double sc[4];
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Create the vector (val,val,val) */
SIRE_ALWAYS_INLINE Vector::Vector(double val)
{
    for (int i=0; i<3; i++)
        sc[i] = val;

    sc[3] = 0.0;
}

/** Create the vector (xpos,ypos,zpos) */
SIRE_ALWAYS_INLINE Vector::Vector(double x, double y, double z)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
    sc[3] = 0.0;
}

/** Destructor */
SIRE_ALWAYS_INLINE Vector::~Vector()
{}

/** Copy assignment operator */
SIRE_ALWAYS_INLINE const Vector& Vector::operator=(const Vector &other)
{
    if (this != &other)
    {
        for (int i=0; i<4; ++i)
            sc[i] = other.sc[i];
    }

    return *this;
}

/** Comparison operator */
SIRE_ALWAYS_INLINE bool Vector::operator==(const Vector &other) const
{
    return &other == this or
           (sc[0] == other.sc[0] and sc[1] == other.sc[1] and
            sc[2] == other.sc[2]);
}

/** Comparison operator */
SIRE_ALWAYS_INLINE bool Vector::operator!=(const Vector &other) const
{
    return &other != this and
           (sc[0] != other.sc[0] or sc[1] != other.sc[1] or
            sc[2] != other.sc[2]);
}

/** Return a raw pointer to the array of coordinates */
SIRE_ALWAYS_INLINE double* Vector::data()
{
    return &(sc[0]);
}

/** Return a raw pointer to the array of coordinates */
SIRE_ALWAYS_INLINE const double* Vector::data() const
{
    return &(sc[0]);
}

/** Return a raw pointer to the array of coordinates */
SIRE_ALWAYS_INLINE const double* Vector::constData() const
{
    return &(sc[0]);
}

/** Return the x component of the vector */
SIRE_ALWAYS_INLINE double Vector::x() const
{
    return sc[0];
}

/** Return the y component of the vector */
SIRE_ALWAYS_INLINE double Vector::y() const
{
    return sc[1];
}

/** Return the z component of the vector */
SIRE_ALWAYS_INLINE double Vector::z() const
{
    return sc[2];
}

/** Return the length of the vector */
SIRE_ALWAYS_INLINE double Vector::length() const
{
    return std::sqrt( pow_2(sc[0]) + pow_2(sc[1]) + pow_2(sc[2]) );
}

/** Return the length^2 of the vector */
SIRE_ALWAYS_INLINE double Vector::length2() const
{
    return pow_2(sc[0]) + pow_2(sc[1]) + pow_2(sc[2]);
}

/** Return the inverse of the length of the vector */
SIRE_ALWAYS_INLINE double Vector::invLength() const
{
    return double(1) / std::sqrt( pow_2(sc[0]) + pow_2(sc[1]) + pow_2(sc[2]) );
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector& Vector::operator+=(const Vector &other)
{
    for (int i=0; i<3; i++)
        sc[i] += other.sc[i];

    return *this;
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector& Vector::operator-=(const Vector &other)
{
    for (int i=0; i<3; i++)
        sc[i] -= other.sc[i];

    return *this;
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector& Vector::operator*=(const double &val)
{
    for (int i=0; i<3; i++)
        sc[i] *= val;

    return *this;
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE Vector Vector::operator-() const
{
    return Vector(-sc[0],-sc[1],-sc[2]);
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector operator+(const Vector &p1, const Vector &p2)
{
    return Vector(p1.sc[0]+p2.sc[0], p1.sc[1]+p2.sc[1], p1.sc[2]+p2.sc[2]);
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector operator-(const Vector &p1, const Vector &p2)
{
    return Vector(p1.sc[0]-p2.sc[0], p1.sc[1]-p2.sc[1], p1.sc[2]-p2.sc[2]);
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector operator*(const Vector &p1, double c)
{
    return Vector(p1.sc[0]*c, p1.sc[1]*c, p1.sc[2]*c);
}

/** Increment, decrement, negate etc. */
SIRE_ALWAYS_INLINE const Vector operator*(double c, const Vector &p1)
{
    return Vector(p1.sc[0]*c, p1.sc[1]*c, p1.sc[2]*c);
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::Vector)
Q_DECLARE_TYPEINFO(SireMaths::Vector, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Vector )

SIRE_END_HEADER

#endif
