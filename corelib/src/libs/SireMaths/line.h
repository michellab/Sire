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

#ifndef SIREMATHS_LINE_H
#define SIREMATHS_LINE_H

#include "vector.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Line;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::Line&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::Line&);

namespace SireMaths
{

/**
This class represents a line in three-dimensional space. (or two points)
 
@author Christopher Woods
*/
class SIREMATHS_EXPORT Line
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const Line&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, Line&);

public:
    Line();
    Line(const Vector &point0, const Vector &point1);
    ~Line();

    static const char* typeName();
    
    const char* what() const
    {
        return Line::typeName();
    }

    QString toString() const;

    double length() const;
    Vector vector() const;

    int count() const;
    
    const Vector& point(int i) const;
    const Vector& operator[](int i) const;
    const Vector& at(int i) const;

private:

    /** The two points that make up the line */
    Vector points[2];
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the number of elements in this line (2) */
SIRE_ALWAYS_INLINE int Line::count() const
{
    return 2;
}

/** Return the i'th point */
SIRE_ALWAYS_INLINE const Vector& Line::point(int i) const
{
    return points[ i%2 ];
}

/** Return the i'th point */
SIRE_ALWAYS_INLINE const Vector& Line::at(int i) const
{
    return point(i);
}

/** Return the i'th point */
SIRE_ALWAYS_INLINE const Vector& Line::operator[](int i) const
{
    return point(i);
}

/** Return the vector that represents this line (goes from point 0 to point 1) */
SIRE_ALWAYS_INLINE Vector Line::vector() const
{
    return (points[1] - points[0]);
}

/** Return the length of the line */
SIRE_ALWAYS_INLINE double Line::length() const
{
    return vector().length();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::Line)
Q_DECLARE_TYPEINFO(SireMaths::Line, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Line )

SIRE_END_HEADER

#endif
