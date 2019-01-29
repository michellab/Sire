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

#ifndef SIREMATHS_PLANE_H
#define SIREMATHS_PLANE_H

#include "vector.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Plane;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::Plane&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::Plane&);

namespace SireMaths
{

/** This class represents an infinite plane, represented by a normal vector 
    perpendicular to the plane, and the distance from the origin 
    along that normal vector.
 
    @author Christopher Woods
*/
class SIREMATHS_EXPORT Plane
{

friend QDataStream& ::operator<<(QDataStream&, const Plane&);
friend QDataStream& ::operator>>(QDataStream&, Plane&);

public:
    Plane();
    Plane(const Vector &normal, const double &distance);
    Plane(const double &a, const double &b, const double &c, const double &d);
    Plane(const Vector &normal, const Vector &contains_point);
    Plane(const Plane &other);
    ~Plane();

    static const char* typeName();
    
    const char* what() const
    {
        return Plane::typeName();
    }

    const Vector& normal() const;
    const double& distanceFromOrigin() const;
    
    double a() const;
    double b() const;
    double c() const;
    double d() const;

    double distance(const Vector &point) const;

private:

    Vector norm;
    double dist;

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the (normalised) normal vector to the plane */
inline const Vector& Plane::normal() const
{
    return norm;
}

/** Return the distance of the plane from the origin along the normal vector */
inline const double& Plane::distanceFromOrigin() const
{
    return dist;
}

/** Return the 'a' component of the plane equation "ax + by + cz + d = 0" */
inline double Plane::a() const
{
    return norm.x();
}

/** Return the 'b' component of the plane equation "ax + by + cz + d = 0" */
inline double Plane::b() const
{
    return norm.y();
}

/** Return the 'c' component of the plane equation "ax + by + cz + d = 0" */
inline double Plane::c() const
{
    return norm.z();
}

/** Return the 'd' component of the plane equation "ax + by + cz + d = 0" */
inline double Plane::d() const
{
    return dist;
}

/** Return the shortest (signed) distance from the plane to the point 'point' */
inline double Plane::distance(const Vector &point) const
{
    return Vector::dot(norm, point) + dist;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::Plane)
Q_DECLARE_TYPEINFO(SireMaths::Plane, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Plane )

SIRE_END_HEADER

#endif
