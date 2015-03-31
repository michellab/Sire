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

#include "SireStream/datastream.h"

#include "line.h"

using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<Line> r_line(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const SireMaths::Line &line)
{
    writeHeader(ds, r_line, 1) << line.points[0] << line.points[1];

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, SireMaths::Line &line)
{
    VersionID v = readHeader(ds, r_line);

    if (v == 1)
    {
        ds >> line.points[0] >> line.points[1];
    }
    else
        throw version_error(v, "1", r_line, CODELOC);

    return ds;
}

/** Construct a zero line */
Line::Line()
{}

/** Construct a line from point0 to point1 */
Line::Line(const Vector &point0, const Vector &point1)
{
    points[0] = point0;
    points[1] = point1;
}

/** Destructor */
Line::~Line()
{}

/** Return a string representation of the line */
QString Line::toString() const
{
    return QObject::tr("Line from %1 to %2, length = %3")
            .arg(point(0).toString(),point(1).toString()).arg(length());
}

const char* Line::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Line>() );
}
