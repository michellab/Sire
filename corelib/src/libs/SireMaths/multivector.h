/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2016  Christopher Woods
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

#ifndef SIREMATHS_MULTIVECTOR_H
#define SIREMATHS_MULTIVECTOR_H

#include "SireMaths/vector.h"
#include "SireMaths/multidouble.h"

#ifdef SIRE_HAS_CPP_11
    #include <functional>
#endif

SIRE_BEGIN_HEADER

namespace SireMaths
{

class MultiVector;
class MultiQuaternion;

SIREMATHS_EXPORT MultiVector operator+(const MultiVector &p1, const MultiVector &p2);
SIREMATHS_EXPORT MultiVector operator-(const MultiVector &p1, const MultiVector &p2);
SIREMATHS_EXPORT MultiVector operator*(const MultiVector &p1, const MultiDouble &c);
SIREMATHS_EXPORT MultiVector operator*(const MultiDouble &c, const MultiVector &p1);
SIREMATHS_EXPORT MultiVector operator/(const MultiVector &p1, const MultiDouble &c);
SIREMATHS_EXPORT MultiQuaternion operator*(const MultiVector &p1, const MultiQuaternion &p2);
SIREMATHS_EXPORT MultiQuaternion operator*(const MultiQuaternion &p1, const MultiVector &p2);

/**
This is a vectorised version of Vector, e.g. x, y, and z 
are held as MultiDouble objects, meaning that this is a 
packed vector of vectors

@author Christopher Woods
*/
class SIREMATHS_EXPORT MultiVector
{
public:
    typedef MultiDouble value_type;

    MultiVector();
    MultiVector(const MultiDouble &value);
    MultiVector(const MultiDouble &x, const MultiDouble &y, const MultiDouble &z);

    MultiVector(const Vector *array, int size);
    MultiVector(const QVector<Vector> &array);
    
    #ifdef SIRE_HAS_CPP_11
        MultiVector(const std::function<double ()> &func);
        MultiVector(const std::function<Vector ()> &func);
    #endif
    
    MultiVector(const MultiVector &other);
    
    ~MultiVector();

    static QVector<MultiVector> fromArray(const Vector *array, int size);
    static QVector<MultiVector> fromArray(const QVector<Vector> &array);

    static const char* typeName();

    const char* what() const
    {
        return MultiVector::typeName();
    }

    static int size();
    static int count();

    MultiDouble* data();
    const MultiDouble* data() const;
    const MultiDouble* constData() const;

    MultiDouble bearing() const;
    MultiDouble bearingXY(const MultiVector &v) const;
    MultiDouble bearingXZ(const MultiVector &v) const;
    MultiDouble bearingYZ(const MultiVector &v) const;

    const MultiDouble& x() const;
    const MultiDouble& y() const;
    const MultiDouble& z() const;

    MultiDouble r() const;
    MultiDouble g() const;
    MultiDouble b() const;

    MultiVector direction() const;
    MultiDouble magnitude() const;

    static MultiDouble angle(const MultiVector &v0, const MultiVector &v1);
    static MultiDouble angle(const MultiVector &v0, const MultiVector &v1,
                             const MultiVector &v2);

    static MultiDouble dihedral(const MultiVector &v0, const MultiVector &v1,
                                const MultiVector &v2, const MultiVector &v3);

    static MultiVector generate(const MultiDouble &dst, const MultiVector &v1,
                                const MultiDouble &ang,
                                const MultiVector &v2, const MultiDouble &dih,
                                const MultiVector &v3);

    MultiVector& operator=(const MultiVector &other);

    bool operator==(const MultiVector &p1) const;
    bool operator!=(const MultiVector &p1) const;
    
    MultiVector& operator+=(const MultiVector &other);
    MultiVector& operator-=(const MultiVector &other);
    MultiVector& operator*=(const MultiDouble &other);
    MultiVector& operator/=(const MultiDouble &other);
    MultiVector operator-() const;

    void set(const MultiDouble &x, const MultiDouble &y, const MultiDouble &z);
    void setX(const MultiDouble &val);
    void setY(const MultiDouble &val);
    void setZ(const MultiDouble &val);
    void setR(const MultiDouble &val);
    void setG(const MultiDouble &val);
    void setB(const MultiDouble &val);

    void set(int i, const Vector &val);
    void quickSet(int i, const Vector &val);

    Vector operator[](int i) const;
    Vector at(int i) const;
    Vector getitem(int i) const;

    MultiDouble manhattanLength() const;
    
    MultiDouble length() const;
    MultiDouble length2() const;
    
    MultiDouble invLength() const;
    MultiDouble invLength2() const;
    
    MultiVector normalise() const;

    QString toString() const;

    static MultiDouble dot(const MultiVector &v0, const MultiVector &v1);
    static MultiVector cross(const MultiVector &v0, const MultiVector &v1);

    void setMax(const MultiVector &other);
    void setMin(const MultiVector &other);

    MultiVector max(const MultiVector &other) const;
    MultiVector min(const MultiVector &other) const;

    static MultiDouble distance2(const MultiVector &v1, const MultiVector &v2);
    static MultiDouble distance(const MultiVector &v1, const MultiVector &v2);
    
    static MultiDouble invDistance(const MultiVector &v1, const MultiVector &v2);
    static MultiDouble invDistance2(const MultiVector &v1, const MultiVector &v2);

    friend MultiVector operator+(const MultiVector &p1, const MultiVector &p2);
    friend MultiVector operator-(const MultiVector &p1, const MultiVector &p2);
    friend MultiVector operator*(const MultiVector &p1, const MultiDouble &c);
    friend MultiVector operator*(const MultiDouble &c, const MultiVector &p1);
    friend MultiVector operator/(const MultiVector &p1, const MultiDouble &c);
    friend MultiQuaternion operator*(const MultiVector &p1, const MultiQuaternion &p2);
    friend MultiQuaternion operator*(const MultiQuaternion &p1, const MultiVector &p2);

    static void swap(MultiVector &v0, int idx0, MultiVector &v1, int idx1);

protected:
    /** The three values, representing the vectorised x, y and z components */
    MultiDouble sc[3];
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Construct an empty vector */
inline MultiVector::MultiVector()
{}

/** Create the vector (val,val,val) */
inline MultiVector::MultiVector(const MultiDouble &val)
{
    for (int i=0; i<3; i++)
        sc[i] = val;
}

/** Create the vector (xpos,ypos,zpos) */
inline MultiVector::MultiVector(const MultiDouble &x, const MultiDouble &y, const MultiDouble &z)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
}

#ifdef SIRE_HAS_CPP_11
    /** Construct a MultiVector from the passed function */
    inline MultiVector::MultiVector(const std::function<double ()> &generator)
    {
        sc[0] = MultiDouble(generator);
        sc[1] = MultiDouble(generator);
        sc[2] = MultiDouble(generator);
    }
    
    /** Construct a MultiVector from the passed function */
    inline MultiVector::MultiVector(const std::function<Vector ()> &generator)
    {
        for (int i=0; i<MultiDouble::count(); ++i)
        {
            Vector v = generator();
            sc[0].set(i, v.x());
            sc[1].set(i, v.y());
            sc[2].set(i, v.z());
        }
    }
#endif

/** Copy assignment operator */
inline MultiVector& MultiVector::operator=(const MultiVector &other)
{
    for (int i=0; i<3; ++i)
        sc[i] = other.sc[i];

    return *this;
}

/** Comparison operator */
inline bool MultiVector::operator==(const MultiVector &other) const
{
    return &other == this or
           (sc[0] == other.sc[0] and sc[1] == other.sc[1] and
            sc[2] == other.sc[2]);
}

/** Comparison operator */
inline bool MultiVector::operator!=(const MultiVector &other) const
{
    return not operator==(other);
}

/** Return a raw pointer to the array of coordinates */
inline MultiDouble* MultiVector::data()
{
    return &(sc[0]);
}

/** Return a raw pointer to the array of coordinates */
inline const MultiDouble* MultiVector::data() const
{
    return &(sc[0]);
}

/** Return a raw pointer to the array of coordinates */
inline const MultiDouble* MultiVector::constData() const
{
    return &(sc[0]);
}

/** Return the x component of the vector */
inline const MultiDouble& MultiVector::x() const
{
    return sc[0];
}

/** Return the y component of the vector */
inline const MultiDouble& MultiVector::y() const
{
    return sc[1];
}

/** Return the z component of the vector */
inline const MultiDouble& MultiVector::z() const
{
    return sc[2];
}

/** Return the length of the vector */
inline MultiDouble MultiVector::length() const
{
    MultiDouble lgth2(sc[0]);
    lgth2 *= sc[0];
    lgth2.multiplyAdd( sc[1], sc[1] );
    lgth2.multiplyAdd( sc[2], sc[2] );
    return lgth2.sqrt();
}

/** Return the length^2 of the vector */
inline MultiDouble MultiVector::length2() const
{
    MultiDouble lgth2(sc[0]);
    lgth2 *= sc[0];
    lgth2.multiplyAdd( sc[1], sc[1] );
    lgth2.multiplyAdd( sc[2], sc[2] );
    return lgth2;
}

/** Return the inverse of the length of the vector */
inline MultiDouble MultiVector::invLength() const
{
    return length().reciprocal();
}

/** Increment, decrement, negate etc. */
inline MultiVector& MultiVector::operator+=(const MultiVector &other)
{
    for (int i=0; i<3; i++)
        sc[i] += other.sc[i];

    return *this;
}

/** Increment, decrement, negate etc. */
inline MultiVector& MultiVector::operator-=(const MultiVector &other)
{
    for (int i=0; i<3; i++)
        sc[i] -= other.sc[i];

    return *this;
}

/** Increment, decrement, negate etc. */
inline MultiVector& MultiVector::operator*=(const MultiDouble &val)
{
    for (int i=0; i<3; i++)
        sc[i] *= val;

    return *this;
}

/** Increment, decrement, negate etc. */
inline MultiVector MultiVector::operator-() const
{
    return MultiVector(-sc[0],-sc[1],-sc[2]);
}

/** Increment, decrement, negate etc. */
inline MultiVector operator+(const MultiVector &p1, const MultiVector &p2)
{
    return MultiVector(p1.sc[0]+p2.sc[0], p1.sc[1]+p2.sc[1], p1.sc[2]+p2.sc[2]);
}

/** Increment, decrement, negate etc. */
inline MultiVector operator-(const MultiVector &p1, const MultiVector &p2)
{
    return MultiVector(p1.sc[0]-p2.sc[0], p1.sc[1]-p2.sc[1], p1.sc[2]-p2.sc[2]);
}

/** Increment, decrement, negate etc. */
inline MultiVector operator*(const MultiVector &p1, const MultiDouble &c)
{
    return MultiVector(p1.sc[0]*c, p1.sc[1]*c, p1.sc[2]*c);
}

/** Increment, decrement, negate etc. */
inline MultiVector operator*(const MultiDouble &c, const MultiVector &p1)
{
    return MultiVector(p1.sc[0]*c, p1.sc[1]*c, p1.sc[2]*c);
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireMaths::MultiVector )

SIRE_END_HEADER

#endif
