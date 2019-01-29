/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMATHS_DISTVECTOR_H
#define SIREMATHS_DISTVECTOR_H

#include "vector.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class DistVector;
}

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::DistVector&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::DistVector&);

namespace SireMaths
{

/** This is a vector that stores the vector as a unit vector giving
    the direction, and a scalar giving the magnitude
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT DistVector : private Vector
{

friend QDataStream& ::operator<<(QDataStream&, const DistVector&);
friend QDataStream& ::operator>>(QDataStream&, DistVector&);

public:
    DistVector();
    
    DistVector(const Vector &vec);
    
    DistVector(const DistVector &other);
    
    ~DistVector();

    static const char* typeName();

    const char* what() const
    {
        return DistVector::typeName();
    }

    double x() const;
    double y() const;
    double z() const;

    double r() const;
    double g() const;
    double b() const;

    const Vector& direction() const;
    double magnitude() const;

    const DistVector& operator=(const DistVector &other);

    bool operator==(const DistVector &p1) const;
    bool operator!=(const DistVector &p1) const;
    
    const DistVector& operator+=(const DistVector &other);
    const DistVector& operator-=(const DistVector &other);
    const DistVector& operator*=(const double &other);
    const DistVector& operator/=(const double &other);
    DistVector operator-() const;

    double operator[](unsigned int i) const;

    double getitem(int i) const{ return this->operator[](i); }

    unsigned int count() const;
    double at(unsigned int i) const;

    double manhattanLength() const;
    
    double length() const;
    double length2() const;
    
    double invLength() const;
    double invLength2() const;
    
    DistVector normalise() const;

    bool isZero() const;

    QString toString() const;

    static DistVector fromString(const QString &str);

    static double dot(const DistVector &v0, const DistVector &v1);
    static DistVector cross(const DistVector &v0, const DistVector &v1);

    void setMax(const DistVector &other);
    void setMin(const DistVector &other);

    DistVector max(const DistVector &other) const;
    DistVector min(const DistVector &other) const;

    Angle bearing() const;
    Angle bearingXY(const DistVector &v) const;
    Angle bearingXZ(const DistVector &v) const;
    Angle bearingYZ(const DistVector &v) const;

    Matrix metricTensor() const;

    static double distance2(const DistVector &v1, const DistVector &v2);
    static double distance(const DistVector &v1, const DistVector &v2);
    
    static double invDistance(const DistVector &v1, const DistVector &v2);
    static double invDistance2(const DistVector &v1, const DistVector &v2);

    static Angle angle(const DistVector &v0, const DistVector &v1);
    static Angle angle(const DistVector &v0, const DistVector &v1, 
                       const DistVector &v2);

    static Angle dihedral(const DistVector &v0, const DistVector &v1,
                          const DistVector &v2, const DistVector &v3);

    static DistVector generate(double dst, const DistVector &v1, const Angle &ang,
                               const DistVector &v2, const Angle &dih, 
                               const DistVector &v3);
};

}

Q_DECLARE_TYPEINFO(SireMaths::DistVector, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(SireMaths::DistVector);

SIRE_EXPOSE_CLASS( SireMaths::DistVector )

SIRE_END_HEADER

#endif
