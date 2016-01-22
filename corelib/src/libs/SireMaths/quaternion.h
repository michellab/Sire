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

#ifndef SIREMATHS_QUATERNION_H
#define SIREMATHS_QUATERNION_H

class QString;

#include "maths.h"

#include "SireUnits/dimensions.h"

#include <QVector>

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Quaternion;
}

class QDataStream;
QDataStream& operator<<(QDataStream&, const SireMaths::Quaternion&);
QDataStream& operator>>(QDataStream&, SireMaths::Quaternion&);

namespace SireMaths
{

class Vector;
class Matrix;

class MultiQuaternion;

const Quaternion operator+(const Quaternion &p1, const Quaternion &p2);
const Quaternion operator-(const Quaternion &p1, const Quaternion &p2);
const Quaternion operator*(const Quaternion &p1, const Quaternion &p2);
const Quaternion operator*(const Quaternion &p1, const Vector &p2);
const Quaternion operator*(const Vector &p1, const Quaternion &p2);

/**
This is a quaternion class that is used to handle 3D rotations and SLERP.
 
@author Christopher Woods
*/
class SIREMATHS_EXPORT Quaternion
{

friend class MultiQuaternion;

friend QDataStream& ::operator<<(QDataStream&, const Quaternion&);
friend QDataStream& ::operator>>(QDataStream&, Quaternion&);

public:
    Quaternion();
    Quaternion(const Quaternion& p1);
    
    Quaternion(SireUnits::Dimension::Angle angle, const Vector &axis);
    Quaternion(const Matrix &m);
    Quaternion(double x, double y, double z, double w);
    
    ~Quaternion();
    
    static const char* typeName();
    
    const char* what() const
    {
        return Quaternion::typeName();
    }

    double x() const;
    double y() const;
    double z() const;
    double w() const;

    bool isIdentity() const;

    Quaternion inverse() const;
    Quaternion conjugate() const;

    double dot(const Quaternion &q) const;

    QString toString() const;
    static Quaternion fromString(const QString &str);

    Matrix toMatrix() const;
    void fromMatrix(const Matrix &m);

    Vector rotate(const Vector &p) const;
    QVector<Vector> rotate(const QVector<Vector> &points) const;

    Quaternion slerp(const Quaternion &q, double lambda) const;

    Quaternion pow(double n) const;

    static Quaternion identity();
    
    void renormalise();

    bool operator==(const Quaternion &p1) const;
    bool operator!=(const Quaternion &p1) const;

    Quaternion& operator=(const Quaternion &p);

    Quaternion& operator+=(const Quaternion &p);
    Quaternion& operator-=(const Quaternion &p);
    Quaternion& operator*=(const Quaternion &p);
    Quaternion& operator*=(const Vector &p);

    friend const Quaternion operator+(const Quaternion &p1, const Quaternion &p2);
    friend const Quaternion operator-(const Quaternion &p1, const Quaternion &p2);
    friend const Quaternion operator*(const Quaternion &p1, const Quaternion &p2);
    friend const Quaternion operator*(const Quaternion &p1, const Vector &p2);
    friend const Quaternion operator*(const Vector &p1, const Quaternion &p2);

private:
    /** The x,y,z,w of the quaternion */
    double sc[4];

};

}

Q_DECLARE_METATYPE(SireMaths::Quaternion)
Q_DECLARE_TYPEINFO(SireMaths::Quaternion, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Quaternion )

SIRE_END_HEADER

#endif
