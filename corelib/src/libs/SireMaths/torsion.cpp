/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "torsion.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireUnits;
using namespace SireMaths;

static RegisterMetaType<Torsion> r_torsion(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const Torsion &torsion)
{
    writeHeader(ds, r_torsion, 1) << torsion.points[0] << torsion.points[1]
                                  << torsion.points[2] << torsion.points[3];

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, Torsion &torsion)
{
    VersionID v = readHeader(ds, r_torsion);

    if (v == 1)
    {
        ds >> torsion.points[0] >> torsion.points[1]
           >> torsion.points[2] >> torsion.points[3];
    }
    else
        throw version_error(v, "1", r_torsion, CODELOC);

    return ds;
}

/** Construct a zero torsion */
Torsion::Torsion()
{}

/** Construct a torsion from the points 0 to 4 */
Torsion::Torsion(const Vector &point0, const Vector &point1,
                 const Vector &point2, const Vector &point3)
{
    points[0] = point0;
    points[1] = point1;
    points[2] = point2;
    points[3] = point3;
}

/** Destructor */
Torsion::~Torsion()
{}

/** Return a string representation of this torsion */
QString Torsion::toString() const
{
    return QObject::tr("Torsion: Angle %1 degrees, length03 = %2")
                  .arg(angle().to(degrees)).arg(vector03().length());
}

/** Return the torsion angle of this torsion (the torsion angle 0-1-2-3 
    around the 1-2 line) */
Angle Torsion::angle() const
{
    return Vector::dihedral(points[0], points[1], points[2], points[3]);
}

/** Return the improper angle of this torsion (the acute angle between the 
    vector 0-1 and the plane formed by 1-2-3) */
Angle Torsion::improperAngle() const
{
    //get the vector perpendicular to the plane
    Vector perp = Vector::cross( points[2] - points[1], points[3] - points[1] );
    
    //get the angle between the perpendicular and the vector
    //from 0-1
    Angle angle = Vector::angle( points[0] - points[1], perp );
    
    //return 90 - angle
    return (90*degrees) - angle;
}

/** Return the line from point 0 to point 3 */
Line Torsion::line03() const
{
    return Line(points[0], points[3]);
}

/** Return the line from point 1 to point 2 */
Line Torsion::line12() const
{
    return Line(points[1], points[2]);
}

/** Return the vector from point 0 to point 3 */
Vector Torsion::vector03() const
{
    return line03().vector();
}

/** Return the vector from point 1 to point 2 */
Vector Torsion::vector12() const
{
    return line12().vector();
}

/** Return the triangle around point 1, i.e. point0-point1-point2 */
Triangle Torsion::triangle1() const
{
    return Triangle(points[0], points[1], points[2]);
}

/** Return the triangle around point 2, i.e. point1-point2-point3 */
Triangle Torsion::triangle2() const
{
    return Triangle(points[1], points[2], points[3]);
}

/** Return the number of points in a torsion (4) */
int Torsion::count() const
{
    return 4;
}

/** Return the point at index 'i' */
const Vector& Torsion::point( int i ) const
{
    return points[ i % 4 ];
}

/** Return the point at index 'i' */
const Vector& Torsion::operator[] ( int i ) const
{
    return this->point(i);
}

/** Return the point at index 'i' */
const Vector& Torsion::at( int i ) const
{
    return this->point(i);
}

const char* Torsion::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Torsion>() );
}
