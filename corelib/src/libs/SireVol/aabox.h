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

#ifndef SIREVOL_AABOX_H
#define SIREVOL_AABOX_H

#include <QVector>

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class AABox;
}

class QDataStream;
SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::AABox&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::AABox&);

namespace SireMaths
{
class Sphere;
}

namespace SireVol
{

using SireMaths::Vector;
using SireMaths::Sphere;

class CoordGroupBase;
class CoordGroupArray;
class CoordGroupArrayArray;

/**
An AABox is an axis-aligned bounding box that is the smallest box that is aligned with the three cartesian axes that completely encases a CoordGroup. It is trivial to obtain the bounding sphere from the AABox. The AABox is used by the distance calculators to quickly determine whether two CoordGroups are within the cutoff radius, and to obtain all CoordGroups that are within particular regions of space.

@author Christopher Woods
*/
class SIREVOL_EXPORT AABox
{

friend QDataStream& ::operator<<(QDataStream&, const AABox&);
friend QDataStream& ::operator>>(QDataStream&, AABox&);

friend class CoordGroupPvt;

public:
    AABox();
    AABox(const Vector &point);
    AABox(const Vector &cent, const Vector &extents);

    AABox(const QVector<Vector> &coordinates);
    AABox(const Vector *coords, int ncoords);

    AABox(const CoordGroupBase &coordgroup);
    AABox(const CoordGroupArray &cgarray);
    AABox(const CoordGroupArrayArray &cgarrays);

    ~AABox();

    static const char* typeName();
    
    const char* what() const
    {
        return AABox::typeName();
    }

    const AABox& operator=(const AABox &other);

    bool operator==(const AABox &other) const;
    bool operator!=(const AABox &other) const;

    AABox& operator+=(const AABox &other);
    AABox& operator+=(const Vector &point);
    AABox& operator+=(const QVector<Vector> &points);

    AABox operator+(const AABox &other) const;
    AABox operator+(const Vector &point) const;
    AABox operator+(const QVector<Vector> &points) const;

    QString toString() const;

    bool isEmpty() const;
    bool isNull() const;

    void add(const AABox &other);
    void add(const Vector &point);
    void add(const QVector<Vector> &points);

    void recalculate(const CoordGroupBase &coordgroup);
    void recalculate(const CoordGroupArray &cgarray);
    void recalculate(const CoordGroupArrayArray &cgarrays);
    void recalculate(const QVector<Vector> &coordinates);

    void translate(const Vector &delta);

    const Vector& center() const;
    const Vector& halfExtents() const;
    Vector maxCoords() const;
    Vector minCoords() const;

    double radius() const;

    Sphere boundingSphere() const;

    bool withinDistance(double dist, const AABox &box) const;
    bool intersects(const AABox &other) const;

    bool contains(const AABox &other) const;
    bool contains(const Vector &point) const;

    static AABox from(const Vector &point);
    static AABox from(const CoordGroupBase &coordgroup);
    static AABox from(const CoordGroupArray &cgarray);
    static AABox from(const CoordGroupArrayArray &cgarrays);
    static AABox from(const QVector<Vector> &coords);
    static AABox from(const Vector &mincoords, const Vector &maxcoords);

protected:
    void recalculate(const Vector *coords, int size);
    void recalculate(const AABox *aaboxes, int size);

private:

    /** The coordinates of the center of this box */
    Vector cent;

    /** The positive half-extents of this box along the x/y/z axes.
        The volume of this box runs from cent-halfextent to cent+halfextent */
    Vector halfextents;

    /** The radius of the smallest sphere that completely contains this box */
    double rad;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Copy operator */
inline const AABox& AABox::operator=(const AABox &box)
{
    cent = box.cent;
    halfextents = box.halfextents;
    rad = box.rad;
    return *this;
}

/** Return the center of the box */
inline const Vector& AABox::center() const
{
    return cent;
}

/** Return the positive half extents of the box */
inline const Vector& AABox::halfExtents() const
{
    return halfextents;
}

/** Return the maximum coordinates of the box */
inline Vector AABox::maxCoords() const
{
    return cent + halfextents;
}

/** Return the minimum coordinates of the box */
inline Vector AABox::minCoords() const
{
    return cent - halfextents;
}

/** Return the radius of the smallest sphere that contains this box
    (the sphere is centered at 'center()', just as the box is) */
inline double AABox::radius() const
{
    return rad;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireVol::AABox)
Q_DECLARE_TYPEINFO(SireVol::AABox, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireVol::AABox )

SIRE_END_HEADER

#endif
