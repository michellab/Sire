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

#include "distvector.h"
#include "matrix.h"

#include "SireUnits/units.h"
#include "SireStream/datastream.h"

#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include <QString>
#include <QRegExp>

#include <cmath>

using namespace SireMaths;
using namespace SireUnits;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<DistVector> r_distvector(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, 
                                         const SireMaths::DistVector &vec)
{
    writeHeader(ds, r_distvector, 1) << vec.sc[0] << vec.sc[1]
                                     << vec.sc[2] << vec.sc[3];

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, 
                                         SireMaths::DistVector &vec)
{
    VersionID v = readHeader(ds, r_distvector);

    if (v == 1)
    {
        ds >> vec.sc[0] >> vec.sc[1] >> vec.sc[2] >> vec.sc[3];
    }
    else
        throw version_error(v, "1", r_distvector, CODELOC);

    return ds;
}

/** Create the zero vector */
DistVector::DistVector() : Vector()
{
    for (int i=0; i<4; i++)
        sc[i] = 0;
}

/** Create from the passed vector */
DistVector::DistVector( const Vector &vec ) : Vector()
{
    double dist = vec.length();
    
    if (SireMaths::isZero(dist))
    {
        for (int i=0; i<4; ++i)
            sc[i] = 0;
    }
    else
    {
        sc[3] = dist;
        dist = 1 / dist;
        
        sc[0] = vec.x() * dist;
        sc[1] = vec.y() * dist;
        sc[2] = vec.z() * dist;
    }
}

/** Construct a DistVector from the QString representation returned by 'toString()' 

    \throw SireError::invalid_arg
*/
DistVector DistVector::fromString(const QString &str)
{
    return DistVector( Vector::fromString(str) );
}

/** Destructor */
DistVector::~DistVector()
{}


/** Copy constructor */
DistVector::DistVector(const DistVector& other)
{
    quickCopy<double>(sc, other.sc, 4);
}

/** Copy assignment operator */
const DistVector& DistVector::operator=(const DistVector &other)
{
    quickCopy<double>(sc, other.sc, 4);

    return *this;
}

/** Comparison operator */
bool DistVector::operator==(const DistVector &other) const
{
    return &other == this or
           (sc[0] == other.sc[0] and sc[1] == other.sc[1] and
            sc[2] == other.sc[2] and sc[3] == other.sc[3]);
}

/** Comparison operator */
bool DistVector::operator!=(const DistVector &other) const
{
    return &other != this and
           (sc[0] != other.sc[0] or sc[1] != other.sc[1] or
            sc[2] != other.sc[2] or sc[3] != other.sc[3]);
}

/** Return the x component of the vector */
double DistVector::x() const
{
    return sc[3] * sc[0];
}

/** Return the y component of the vector */
double DistVector::y() const
{
    return sc[3] * sc[1];
}

/** Return the z component of the vector */
double DistVector::z() const
{
    return sc[3] * sc[2];
}

/** Return the length of the vector */
double DistVector::length() const
{
    return sc[3];
}

/** Return the length^2 of the vector */
double DistVector::length2() const
{
    return pow_2(sc[3]);
}

/** Return the inverse of the length of the vector */
double DistVector::invLength() const
{
    return 1 / sc[3];
}

/** Return the inverse length squared */
double DistVector::invLength2() const
{
    return 1 / pow_2(sc[3]);
}

/** Return the distance squared between two vectors */
double DistVector::distance2(const DistVector &v1, const DistVector &v2)
{
    return Vector::distance2( v1, v2 );
}

/** Return the distance between two vectors */
double DistVector::distance(const DistVector &v1, const DistVector &v2)
{
    return Vector::distance( v1, v2 );
}

/** Return the 1 / distance between two vectors */
double DistVector::invDistance(const DistVector &v1, const DistVector &v2)
{
    return Vector::invDistance( v1, v2 );
}

/** Return 1 / distance2 between two vectors */
double DistVector::invDistance2(const DistVector &v1, const DistVector &v2)
{
    return Vector::invDistance2( v1, v2 );
}

/** Access the elements of the array via an index operator */
double DistVector::operator[](unsigned int i) const
{
    return sc[3] * sc[i%3];
}

/** Return the size of the Vector (always 3 - unless you disagree
    with me that we should be living in a 3-dimensional space!) */
unsigned int DistVector::count() const
{
    return 3;
}

/** Access elements by index */
double DistVector::at(unsigned int i) const
{
    return this->operator[](i);
}

/** Return whether or not this is a zero length vector */
bool DistVector::isZero() const
{
    return SireMaths::isZero(sc[3]) or Vector(*this).isZero();
}

/** Return a normalised form of the vector */
DistVector DistVector::normalise() const
{
    DistVector ret(*this);
    ret.sc[3] = 1;
    return ret;
}

/** Return the dot product of v0 and v1 */
double DistVector::dot(const DistVector &v0, const DistVector &v1)
{
    return (v0.sc[0]*v1.sc[0] + v0.sc[1]*v1.sc[1] + v0.sc[2]*v1.sc[2]);
}

/** Return a QString representation of the vector */
QString DistVector::toString() const
{
    return QObject::tr("( %1, %2, %3 )").arg(x()).arg(y()).arg(z());
}

/** Return the components via rgb (limited between 0 and 1) */
double DistVector::r() const
{
    return std::max(0.0, std::min(1.0,x()));
}

/** Return the components via rgb (limited between 0 and 1) */
double DistVector::g() const
{
    return std::max(0.0, std::min(1.0,y()));
}

/** Return the components via rgb (limited between 0 and 1) */
double DistVector::b() const
{
    return std::max(0.0, std::min(1.0,z()));
}

/** Return the direction of this vector */
const Vector& DistVector::direction() const
{
    return *this;
}

/** Return the magnitude of this vector */
double DistVector::magnitude() const
{
    return sc[3];
} 

/** Return the bearing of this vector against (0,1,0) (north) on the xy plane */
Angle DistVector::bearing() const
{
    return this->direction().bearing();
}

/** Return the bearing of this vector against 'v' on the xy plane */
Angle DistVector::bearingXY(const DistVector &v) const
{
    return this->direction().bearingXY( v.direction() );
}

/** Return the bearing of this vector against 'v' on the xz plane */
Angle DistVector::bearingXZ(const DistVector &v) const
{
    return this->direction().bearingXZ( v.direction() );
}

/** Return the bearing of this vector against 'v' on the yz plane */
Angle DistVector::bearingYZ(const DistVector &v) const
{
    return this->direction().bearingYZ( v.direction() );
}

/** Set this Vector so that it has the maximum x/y/z components out of
    this and 'other' (e.g. this->x = max(this->x(),other.x() etc.) */
void DistVector::setMax(const DistVector &other)
{
    this->operator=( Vector(*this).max( Vector(other) ) );
}

/** Set this Vector so that it has the minimum x/y/z components */
void DistVector::setMin(const DistVector &other)
{
    this->operator=( Vector(*this).min( Vector(other) ) );
}

/** Return a vector that has the maximum x/y/z components out of this
    and 'other' */
DistVector DistVector::max(const DistVector &other) const
{
    DistVector ret(*this);
    ret.setMax(other);
    return ret;
}

/** Return a vector that has the minimum components */
DistVector DistVector::min(const DistVector &other) const
{
    DistVector ret(*this);
    ret.setMin(other);
    return ret;
}

/** Return the angle between vectors 'v0' and 'v1' - this is the smallest
    angle, and will always lie between 0 and 180 degrees */
Angle DistVector::angle(const DistVector &v0, const DistVector &v1)
{
    return Vector::angle( v0.direction(), v1.direction() );
}

/** Return the angle between v0-v1-v2 (treating the vectors as points in space) */
Angle DistVector::angle(const DistVector &v0, const DistVector &v1, const DistVector &v2)
{
    return Vector::angle( v0, v1, v2 );
}

/** Return the dihedral angle between v0-v1-v2-v3 (treating the vectors as points) */
Angle DistVector::dihedral(const DistVector &v0, const DistVector &v1, const DistVector &v2, const DistVector &v3)
{
    return Vector::dihedral( v0, v1, v2, v3 );
}
    
/** Generate a vector, v0, that has distance 'dst' v0-v1, angle 'ang' v0-v1-v2,
    and dihedral 'dih' v0-v1-v2-v3 */
DistVector DistVector::generate(double dst, const DistVector &v1, const Angle &ang, 
                                const DistVector &v2,
                                const Angle &dih, const DistVector &v3)
{
    return DistVector( Vector::generate(dst, v1, ang, v2, dih, v3) );
}

 /** Return the cross product of v0 and v1 */
DistVector DistVector::cross(const DistVector &v0, const DistVector &v1)
{
    return Vector::cross( v0.direction(), v1.direction() );
}

/** Return the manhattan length of the vector */
double DistVector::manhattanLength() const
{
    return Vector(*this).manhattanLength();
}

/** Return the metric tensor of a vector, i.e.
          
    | y*y + z*z,    -x*y    -x*z      |
    |    -y*x,   x*x + z*z  -y*z      |
    |    -z*x       -z*y    x*x + y*y |

*/
Matrix DistVector::metricTensor() const
{
    return Vector(*this).metricTensor();
}

/** Increment, decrement, negate etc. */
const DistVector& DistVector::operator+=(const DistVector &other)
{
    return this->operator=( Vector(*this) + Vector(other) );
}

/** Increment, decrement, negate etc. */
const DistVector& DistVector::operator-=(const DistVector &other)
{
    return this->operator=( Vector(*this) - Vector(other) );
}

/** Increment, decrement, negate etc. */
const DistVector& DistVector::operator*=(const double &val)
{
    sc[3] *= val;
    return *this;
}

/** Increment, decrement, negate etc. */
const DistVector& DistVector::operator/=(const double &val)
{
    if ( SireMaths::isZero(val) )
        throw SireMaths::math_error(QObject::tr(
            "Cannot divide a vector by zero! %1 / 0 is a error!").arg(this->toString()),CODELOC);
    
    sc[3] /= val;
    return *this;
}

/** Increment, decrement, negate etc. */
DistVector DistVector::operator-() const
{
    DistVector ret(*this);
    ret.sc[3] = -sc[3];
    return ret;
}

const char* DistVector::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DistVector>() );
}
