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

#include "aabox.h"
#include "coordgroup.h"

#include "SireMaths/sphere.h"

#include <QDebug>

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireVol;

static const RegisterMetaType<AABox> r_aabox(NO_ROOT);

/** Serialise an AABox to a binary datastream */
QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const AABox &aabox)
{
    writeHeader(ds, r_aabox, 1) << aabox.cent << aabox.halfextents << aabox.rad;

    return ds;
}

/** Deserialise an AABox from a binary datastream */
QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, AABox &aabox)
{
    VersionID v = readHeader(ds, r_aabox);

    if (v == 1)
    {
        ds >> aabox.cent >> aabox.halfextents >> aabox.rad;
    }
    else
        throw version_error(v, "1", r_aabox, CODELOC);

    return ds;
}

/** Construct an empty AABox */
AABox::AABox() : cent(), halfextents(), rad(0)
{}

/** Construct an AABox that completely encloses the point 'point' */
AABox::AABox(const Vector &point) : cent(point), halfextents(), rad(0)
{}

/** Construct an AABox with center at 'cent', and half-extents 'extents' */
AABox::AABox(const Vector &c, const Vector &extents) : cent(c), halfextents(extents)
{
    rad = halfextents.length();
}

/** Construct an AABox that completely encases the CoordGroup 'coordgroup' */
AABox::AABox(const CoordGroupBase &coordgroup)
{
    recalculate(coordgroup);
}

/** Construct an AABox that completely encases all of the points in all of the
    CoordGroups in 'cgarray' */
AABox::AABox(const CoordGroupArray &cgarray)
{
    recalculate(cgarray.constAABoxData(), cgarray.nCoordGroups());
}

/** Construct an AABox that completely encases all of the points in all of the
    CoordGroups in all of the arrays in 'cgarrays' */
AABox::AABox(const CoordGroupArrayArray &cgarrays)
{
    recalculate(cgarrays.constAABoxData(), cgarrays.nCoordGroups());
}

/** Construct an AABox that completely encases the points  in 'coordinates' */
AABox::AABox(const QVector<Vector> &coordinates)
{
    recalculate(coordinates);
}

/** Construct an AABox that completely encases the points in 'coords' */
AABox::AABox(const Vector *coords, int ncoords)
{
    recalculate(coords, ncoords);
}

/** Destructor */
AABox::~AABox()
{}

/** Comparison operator */
bool AABox::operator==(const AABox &other) const
{
    return this == &other or
          (rad == other.rad and cent == other.cent and halfextents == other.halfextents);
}

/** Comparison operator */
bool AABox::operator!=(const AABox &other) const
{
    return this != &other and
          (rad != other.rad or cent != other.cent or halfextents != other.halfextents);
}

/** Return if the AABox is null */
bool AABox::isNull() const
{
    return cent.isZero() and rad == 0;
}

/** Return if the AABox is empty */
bool AABox::isEmpty() const
{
    return rad == 0;
}

/** Return an AABox constructed to contain the coordinates in 'coordinates' */
AABox AABox::from(const QVector<Vector> &coordinates)
{
    return AABox(coordinates);
}

/** Return an AABox constructed to contain the coordinates of 'point' */
AABox AABox::from(const Vector &point)
{
    return AABox(point);
}

/** Return an AABox constructed to contain the coordinates in 'coordinates' */
AABox AABox::from(const CoordGroupBase &coordinates)
{
    return AABox(coordinates);
}

/** Return an AABox constructed to contain all of the CoordGroups in 'cgarray' */
AABox AABox::from(const CoordGroupArray &cgarray)
{
    return AABox(cgarray);
}

/** Return an AABox constructed to contain all of the CoordGroups in the 
    arrays in 'cgarrays' */
AABox AABox::from(const CoordGroupArrayArray &cgarrays)
{
    return AABox(cgarrays);
}

/** Return a string representation of this AABox */
QString AABox::toString() const
{
    return QObject::tr( "AABox( min=%1, max=%2 )" )
                .arg(this->minCoords().toString(), this->maxCoords().toString() );
}

/** Internal function used to recalculate the AABox from the passed 
    array of AABoxes */
void AABox::recalculate(const AABox *aaboxes, int sz)
{
    if (sz == 1)
    {
        this->operator=(aaboxes[0]);
    }
    else if (sz > 1)
    {
        //set the initial max and min coords from the first coordinate in the group
        Vector maxcoords( aaboxes[0].maxCoords() );
        Vector mincoords( aaboxes[0].minCoords() );

        //loop through all of the remaining coordinates in the group
        for (int i=1; i < sz; ++i)
        {
            //calculate the maximum and minimum coordinates
            const AABox &aabox = aaboxes[i];
            maxcoords.setMax( aabox.maxCoords() );
            mincoords.setMin( aabox.minCoords() );
        }

        //now calculate the center as half the maximum and minimum coordinates
        cent = 0.5 * (maxcoords + mincoords);

        //the positive half-extent is the difference between the maximum
        //coordinates and the center
        halfextents = maxcoords - cent;

        //the radius is the length of 'halfextents'
        rad = halfextents.length();
    }
    else
    {
        cent = Vector(0);
        halfextents = Vector(0);
        sz = 0;
    }
}

/** Internal function used to recalculate the AABox from the coordinates in the
    array 'coords' (which has size 'sz') */
void AABox::recalculate(const Vector *coords, int sz)
{
    if (sz > 0)
    {
        //set the initial max and min coords from the first coordinate in the group
        Vector maxcoords( coords[0] );
        Vector mincoords( maxcoords );

        //loop through all of the remaining coordinates in the group
        for (int i=1; i < sz; ++i)
        {
            //calculate the maximum and minimum coordinates
            const Vector &coord = coords[i];
            maxcoords.setMax( coord );
            mincoords.setMin( coord );
        }

        //now calculate the center as half the maximum and minimum coordinates
        cent = 0.5 * (maxcoords + mincoords);

        //the positive half-extent is the difference between the maximum
        //coordinates and the center
        halfextents = maxcoords - cent;

        //the radius is the length of 'halfextents'
        rad = halfextents.length();
    }
    else
    {
        cent = Vector(0);
        halfextents = Vector(0);
        sz = 0;
    }
}

/** Recalculate the AABox so that it completely encloses the CoordGroup 'coordgroup' */
void AABox::recalculate(const CoordGroupBase &coordgroup)
{
    this->recalculate( coordgroup.constData(), coordgroup.size() );
}

/** Recalculate the AABox so that it completely encloses the CoordGroups
    in the array 'cgarray' */
void AABox::recalculate(const CoordGroupArray &cgarray)
{
    this->recalculate( cgarray.constAABoxData(), cgarray.nCoordGroups() );
}

/** Recalculate the AABox so that it completely encloses the CoordGroups
    in the arrays 'cgarrays' */
void AABox::recalculate(const CoordGroupArrayArray &cgarrays)
{
    this->recalculate( cgarrays.constAABoxData(), cgarrays.nCoordGroups() );
}

/** Recalculate the AABox so that it completely encloses the 'coordinates' */
void AABox::recalculate(const QVector<Vector> &coordinates)
{
    this->recalculate( coordinates.constData(), coordinates.size() );
}

/** Return whether or not this box is within 'dist' of box 'box'.
    (using infinite cartesian axes) */
bool AABox::withinDistance(double dist, const AABox &box) const
{
   //look at the components of the distance along the x, y and z axes
    double dx = std::abs(cent.x() - box.cent.x()) - halfextents.x() - box.halfextents.x();
    double dy = std::abs(cent.y() - box.cent.y()) - halfextents.y() - box.halfextents.y();
    double dz = std::abs(cent.z() - box.cent.z()) - halfextents.z() - box.halfextents.z();

    dx = SIRE_MAX(dx,0.0);
    dy = SIRE_MAX(dy,0.0);
    dz = SIRE_MAX(dz,0.0);

    return dx*dx + dy*dy + dz*dz <= dist*dist;
}

/** Return whether this box intersects with 'box' */
bool AABox::intersects(const AABox &box) const
{
    //look at the components of the distance along the x, y and z axes
    double dx = std::abs(cent.x() - box.cent.x()) - halfextents.x() - box.halfextents.x();
    double dy = std::abs(cent.y() - box.cent.y()) - halfextents.y() - box.halfextents.y();
    double dz = std::abs(cent.z() - box.cent.z()) - halfextents.z() - box.halfextents.z();

    return dx <= 0.0 and dy <= 0.0 and dz <= 0.0;
}

/** Return whether or not this box contains 'other' */
bool AABox::contains(const AABox &other) const
{
    const Vector mindelta = this->minCoords() - other.minCoords();
    
    if (mindelta.x() > 0 or mindelta.y() > 0 or mindelta.z() > 0)
        return false;

    else
    {
        const Vector maxdelta = this->maxCoords() - other.maxCoords();
        return (maxdelta.x() > 0 and maxdelta.y() > 0 and maxdelta.z() > 0);
    }
}

/** Return whether or not this box contains the point 'point' */
bool AABox::contains(const Vector &point) const
{
    const Vector mindelta = this->minCoords() - point;
    
    if (mindelta.x() > 0 or mindelta.y() > 0 or mindelta.z() > 0)
        return false;
    
    else
    {
        const Vector maxdelta = this->maxCoords() - point;
        return (maxdelta.x() > 0 and maxdelta.y() > 0 and maxdelta.z() > 0);
    }
}

/** Translate this AABox by 'delta' */
void AABox::translate(const Vector &delta)
{
    cent += delta;
}

/** Add another AABox to this one - this forms the union of both of the
    boxes. */
AABox& AABox::operator+=(const AABox &other)
{
    Vector mincoords = this->minCoords();
    Vector maxcoords = this->maxCoords();

    Vector old_mincoords = mincoords;
    Vector old_maxcoords = maxcoords;

    mincoords.setMin( other.minCoords() );
    maxcoords.setMax( other.maxCoords() );

    if (mincoords != old_mincoords or
        maxcoords != old_maxcoords)
    {
        cent = 0.5 * (maxcoords + mincoords);

        halfextents = maxcoords - cent;

        rad = halfextents.length();
    }

    return *this;
}

/** Add another AABox to this one - this forms the union of both of the
    boxes. */
AABox AABox::operator+(const AABox &other) const
{
    AABox ret = *this;
    ret += other;
    return ret;
}

/** Add another AABox to this one - this forms the union of both of the
    boxes. */
void AABox::add(const AABox &other)
{
    *this += other;
}

/** Add a point to this box */
AABox& AABox::operator+=(const Vector &point)
{
    Vector mincoords = this->minCoords();
    Vector maxcoords = this->maxCoords();

    Vector old_mincoords = mincoords;
    Vector old_maxcoords = maxcoords;

    mincoords.setMin( point );
    maxcoords.setMax( point );

    if (old_mincoords != mincoords or
        old_maxcoords != maxcoords)
    {
        cent = 0.5 * (maxcoords + mincoords);

        halfextents = maxcoords - cent;

        rad = halfextents.length();
    }

    return *this;
}

/** Add a point to this box */
AABox AABox::operator+(const Vector &point) const
{
    AABox ret = *this;
    ret += point;
    return ret;
}

/** Add a point to this box */
void AABox::add(const Vector &point)
{
    *this += point;
}

/** Add lots of points to this box */
AABox& AABox::operator+=(const QVector<Vector> &points)
{
    if (points.isEmpty())
        return *this;

    Vector mincoords = this->minCoords();
    Vector maxcoords = this->maxCoords();

    Vector old_mincoords = mincoords;
    Vector old_maxcoords = maxcoords;

    int npoints = points.count();
    const Vector *points_array = points.constData();

    for (int i=0; i<npoints; ++i)
    {
        const Vector &point = points_array[i];
        mincoords.setMin(point);
        maxcoords.setMax(point);
    }

    if (mincoords != old_mincoords or
        maxcoords != old_maxcoords)
    {
        cent = 0.5 * (maxcoords + mincoords);

        halfextents = maxcoords - cent;

        rad = halfextents.length();
    }

    return *this;
}

/** Construct a new AABox from the passed minimum and maximum coordinates */
AABox AABox::from(const Vector &mincoords, const Vector &maxcoords)
{
    Vector min = mincoords;
    min.setMin(maxcoords);
    
    Vector max = maxcoords;
    max.setMax(mincoords);
    
    Vector halfextents = 0.5 * (max - min);
    Vector cent = min + halfextents;
    
    return AABox(cent, halfextents);
}

/** Add lots of points to this box */
AABox AABox::operator+(const QVector<Vector> &points) const
{
    AABox ret = *this;
    ret += points;
    return ret;
}

/** Add lots of points to this box */
void AABox::add(const QVector<Vector> &points)
{
    *this += points;
}

/** Return the sphere that just contains this AABox */
Sphere AABox::boundingSphere() const
{
    return Sphere( this->center(), this->radius() );
}

const char* AABox::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AABox>() );
}
