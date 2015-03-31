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

#include <QDataStream>

#include "triangle.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireUnits;
using namespace SireMaths;

static const RegisterMetaType<Triangle> r_triangle(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Triangle &triangle)
{
    writeHeader(ds, r_triangle, 1)
          << triangle.points[0] << triangle.points[1] << triangle.points[2];

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Triangle &triangle)
{
    VersionID v = readHeader(ds, r_triangle);

    if (v == 1)
    {
        ds >> triangle.points[0] >> triangle.points[1] >> triangle.points[2];
    }
    else
        throw version_error(v, "1", r_triangle, CODELOC);

    return ds;
}

/** Create a zero triangle */
Triangle::Triangle()
{}

/** Create a triangle from point 0 -> 1 -> 2 */
Triangle::Triangle(const Vector &point0, const Vector &point1, const Vector &point2)
{
    points[0] = point0;
    points[1] = point1;
    points[2] = point2;
}

/** Copy constructor */
Triangle::Triangle(const Triangle &other)
{
    for (int i=0; i<3; ++i)
        points[i] = other.points[i];
}

/** Destructor */
Triangle::~Triangle()
{}

/** Return a string representation of the triangle */
QString Triangle::toString() const
{
    return QObject::tr("Triangle: Angles %1 degs, %2 degs, %3 degs")
                .arg(angle0().to(degrees))
                .arg(angle1().to(degrees))
                .arg(angle2().to(degrees));
}

const char* Triangle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Triangle>() );
}
