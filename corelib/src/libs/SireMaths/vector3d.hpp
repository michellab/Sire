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

#ifndef SIREMATHS_VECTOR3D_HPP
#define SIREMATHS_VECTOR3D_HPP

#include "tostring.h"
#include "vector.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
template<class T>
class Vector3D;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireMaths::Vector3D<T>&);

template<class T>
QDataStream& operator>>(QDataStream&, SireMaths::Vector3D<T>&);

namespace SireMaths
{

/** This is a simple small class that holds 3D value
    (e.g. 3D velocity, 3D forces)

    @author Christopher Woods
*/
template<class T>
class Vector3D
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<<>(QDataStream&, const Vector3D<T>&);
friend SIREMATHS_EXPORT QDataStream& ::operator>><>(QDataStream&, Vector3D<T>&);

public:
    typedef T value_type;

    Vector3D();
    Vector3D(const T &val);
    Vector3D(const T &x, const T &y, const T &z);
    Vector3D(const Vector &v);

    Vector3D(const Vector3D<T> &other);

    ~Vector3D();

    Vector3D<T>& operator=(const Vector3D<T> &other);

    bool operator==(const Vector3D<T> &other) const;
    bool operator!=(const Vector3D<T> &other) const;

    const T& operator[](int i) const;

    Vector3D<T>* clone() const;

    Vector3D<T>& operator+=(const Vector3D<T> &other);
    Vector3D<T>& operator-=(const Vector3D<T> &other);

    Vector3D<T>& operator*=(double scale);
    Vector3D<T>& operator/=(double scale);

    Vector3D<T> operator-() const;
    Vector3D<T> operator+(const Vector3D<T> &other) const;
    Vector3D<T> operator-(const Vector3D<T> &other) const;
    Vector3D<T> operator*(double scale) const;
    Vector3D<T> operator/(double scale) const;

    static const char* typeName();

    const char* what() const;

    const T& x() const;
    const T& y() const;
    const T& z() const;

    void set(int i, const T &value);

    int count() const;

    const T& at(int i) const;

    Vector value() const;

    QString toString() const;

private:
    /** The three values */
    T sc[3];
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor - this default-constructs the values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>::Vector3D()
{}

/** Construct, initialising to 'val' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>::Vector3D(const T &val)
{
    sc[0] = val;
    sc[1] = val;
    sc[2] = val;
}

/** Construct from the passed values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>::Vector3D(const T &x, const T &y, const T &z)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
}

/** Construct from the passed vector - this assumes that
    the vector is already in internal units */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>::Vector3D(const Vector &v)
{
    sc[0] = T(v.x());
    sc[1] = T(v.y());
    sc[2] = T(v.z());
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>::Vector3D(const Vector3D<T> &other)
{
    sc[0] = other.sc[0];
    sc[1] = other.sc[1];
    sc[2] = other.sc[2];
}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>::~Vector3D()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>& Vector3D<T>::operator=(const Vector3D<T> &other)
{
    if (this != &other)
    {
        sc[0] = other.sc[0];
        sc[1] = other.sc[1];
        sc[2] = other.sc[2];
    }

    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Vector3D<T>::operator==(const Vector3D<T> &other) const
{
    return (this == &other) or
           (sc[0] == other.sc[0] and sc[1] == other.sc[1] and sc[2] == other.sc[2]);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Vector3D<T>::operator!=(const Vector3D<T> &other) const
{
    return (this != &other) and
           (sc[0] != other.sc[0] or sc[1] != other.sc[1] or sc[2] != other.sc[2]);
}

/** Return the ith component */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Vector3D<T>::operator[](int i) const
{
    return sc[ i % 3 ];
}

/** Return a clone of this vector */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>* Vector3D<T>::clone() const
{
    return new Vector3D<T>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* Vector3D<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Vector3D<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* Vector3D<T>::what() const
{
    return Vector3D<T>::typeName();
}

/** Return the x component */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Vector3D<T>::x() const
{
    return sc[0];
}

/** Return the y component */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Vector3D<T>::y() const
{
    return sc[1];
}

/** Return the z component */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Vector3D<T>::z() const
{
    return sc[2];
}

/** Return this vector in internal units */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector Vector3D<T>::value() const
{
    return Vector( sc[0].value(), sc[1].value(), sc[2].value() );
}

/** Set the ith component equal to 'value' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Vector3D<T>::set(int i, const T &value)
{
    sc[ i % 3 ] = value;
}

/** Return the number of components */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int Vector3D<T>::count() const
{
    return 3;
}

/** Return the ith component */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Vector3D<T>::at(int i) const
{
    return sc[ i % 3 ];
}

/** Return a string representation of this Vector3D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString Vector3D<T>::toString() const
{
    return QString( "( %1, %2, %3 )" )
                .arg( Sire::toString(sc[0]),
                      Sire::toString(sc[1]),
                      Sire::toString(sc[2]) );
}

/** Addition operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>& Vector3D<T>::operator+=(const Vector3D<T> &other)
{
    sc[0] += other.sc[0];
    sc[1] += other.sc[1];
    sc[2] += other.sc[2];
    return *this;
}

/** Subtraction operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>& Vector3D<T>::operator-=(const Vector3D<T> &other)
{
    sc[0] -= other.sc[0];
    sc[1] -= other.sc[1];
    sc[2] -= other.sc[2];
    return *this;
}

/** Multiplication operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>& Vector3D<T>::operator*=(double scale)
{
    sc[0] *= scale;
    sc[1] *= scale;
    sc[2] *= scale;
    return *this;
}

/** Division operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T>& Vector3D<T>::operator/=(double scale)
{
    sc[0] /= scale;
    sc[1] /= scale;
    sc[2] /= scale;
    return *this;
}

/** Negation operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T> Vector3D<T>::operator-() const
{
    return Vector3D<T>(-sc[0], -sc[1], -sc[2]);
}

/** Addition operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T> Vector3D<T>::operator+(const Vector3D<T> &other) const
{
    return Vector3D<T>(sc[0]+other.sc[0],
                       sc[1]+other.sc[1],
                       sc[2]+other.sc[2]);
}

/** Subtraction operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T> Vector3D<T>::operator-(const Vector3D<T> &other) const
{
    return Vector3D<T>(sc[0]-other.sc[0],
                       sc[1]-other.sc[1],
                       sc[2]-other.sc[2]);
}

/** Multiplication operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T> Vector3D<T>::operator*(double scale) const
{
    return Vector3D<T>(scale*sc[0],
                       scale*sc[1],
                       scale*sc[2]);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T> operator*(double scale, const Vector3D<T> &v)
{
    return v * scale;
}

/** Division operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Vector3D<T> Vector3D<T>::operator/(double scale) const
{
    return Vector3D<T>(sc[0]/scale,
                       sc[1]/scale,
                       sc[2]/scale);
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMaths::Vector3D<T> &vec)
{
    ds << vec.sc[0] << vec.sc[1] << vec.sc[2];
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMaths::Vector3D<T> &vec)
{
    ds >> vec.sc[0] >> vec.sc[1] >> vec.sc[2];
    return ds;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
