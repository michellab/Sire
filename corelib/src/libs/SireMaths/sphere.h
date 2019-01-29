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

#ifndef SIREMATHS_SPHERE_H
#define SIREMATHS_SPHERE_H

#include "vector.h"

#include <QVector>

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Sphere;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::Sphere&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::Sphere&);

namespace SireMaths
{

/**
This class is a mathematical representation of a sphere.
 
@author Christopher Woods
*/
class SIREMATHS_EXPORT Sphere
{

friend QDataStream& ::operator<<(QDataStream&, const Sphere&);
friend QDataStream& ::operator>>(QDataStream&, Sphere&);

public:
    Sphere();
    Sphere(const double &radius);
    Sphere(const Vector &position, const double &radius);
    Sphere(const Sphere &other);
    
    ~Sphere();

    static const char* typeName();
    
    const char* what() const
    {
        return Sphere::typeName();
    }

    bool operator==(const Sphere &other) const;
    bool operator!=(const Sphere &other) const;

    QString toString() const;

    const Vector& position() const;
    const Vector& center() const;
    double radius() const;

    Sphere translate(const Vector &delta) const;

    double volume() const;
    double surfaceArea() const;

    void setPosition(const Vector &position);
    void setCenter(const Vector &center);
    void setRadius(double radius);

    bool intersects(const Sphere &other) const;
    
    bool contains(const Vector &point) const;
    bool contains(const Sphere &other) const;
    
    double intersectionVolume(const Sphere &other) const;
    double intersectionVolume(const Sphere &other0, const Sphere &other1) const;
    
    static double combinedVolume(const QVector<Sphere> &spheres);

    static double combinedVolumeMC(const QVector<Sphere> &spheres,
                                   qint64 nsamples=-1);

private:

    /** The location of the center of the sphere */
    Vector _center;

    /** The radius of the sphere */
    double _radius;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the position of the center of the sphere */
inline const Vector& Sphere::position() const
{
    return _center;
}

/** Return the position of the center of the sphere */
inline const Vector& Sphere::center() const
{
    return _center;
}

/** Return the radius of the sphere */
inline double Sphere::radius() const
{
    return _radius;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::Sphere)
Q_DECLARE_TYPEINFO(SireMaths::Sphere, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Sphere )

SIRE_END_HEADER

#endif
