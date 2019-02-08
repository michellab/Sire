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

#ifndef SIREVOL_PERIODICBOX_H
#define SIREVOL_PERIODICBOX_H

#include "cartesian.h"

#include "SireMaths/vector.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireVol
{
class PeriodicBox;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::PeriodicBox&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::PeriodicBox&);

namespace SireVol
{

using SireMaths::Vector;

/**
A PeriodicBox is a volume  that represents standard periodic boundary conditions
(a 3D box replicated to infinity along all three dimensions).

@author Christopher Woods
*/
class SIREVOL_EXPORT PeriodicBox 
        : public SireBase::ConcreteProperty<PeriodicBox,Cartesian>
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const PeriodicBox&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, PeriodicBox&);

public:
    PeriodicBox();
    PeriodicBox(const Vector &extents);
    PeriodicBox(const Vector &min, const Vector &max);

    PeriodicBox(const PeriodicBox &other);

    ~PeriodicBox();

    PeriodicBox& operator=(const PeriodicBox &other);
    
    bool operator==(const PeriodicBox &other) const;
    bool operator!=(const PeriodicBox &other) const;

    bool isPeriodic() const;
    bool isCartesian() const;

    QString toString() const;

    SireUnits::Dimension::Volume volume() const;
    SpacePtr setVolume(SireUnits::Dimension::Volume volume) const;

    void setDimensions(const Vector &dimensions);
    void setDimensions(const Vector &mincoords, const Vector &maxcoords);
    
    const Vector& dimensions() const;

    Vector minCoords(const Vector &center = Vector(0)) const;
    Vector maxCoords(const Vector &center = Vector(0)) const;

    static const char* typeName();

    double calcDist(const Vector &point0, const Vector &point1) const;
    double calcDist2(const Vector &point0, const Vector &point1) const;

    double calcDist(const CoordGroup &group1, const CoordGroup &group2,
                    DistMatrix &distmat) const;

    double calcDist(const CoordGroup &group, const Vector &point,
                    DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group, const Vector &point,
                     DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group1, const CoordGroup &group2,
                     DistMatrix &distmat) const;

    double calcInvDist(const CoordGroup &group1, const CoordGroup &group2,
                       DistMatrix &distmat) const;

    double calcInvDist2(const CoordGroup &group1, const CoordGroup &group2,
                        DistMatrix &distmat) const;

    DistVector calcDistVector(const Vector &point0, const Vector &point1) const;
    
    double calcDistVectors(const CoordGroup &group1, const CoordGroup &group2,
                           DistVectorMatrix &distmat) const;

    double calcDistVectors(const CoordGroup &group, const Vector &point,
                           DistVectorMatrix &distmat) const;

    SireUnits::Dimension::Angle calcAngle(const Vector &point0,
                                          const Vector &point1,
                                          const Vector &point2) const;

    SireUnits::Dimension::Angle calcDihedral(const Vector &point0,
                                             const Vector &point1,
                                             const Vector &point2,
                                             const Vector &point3) const;

    bool beyond(double dist, const AABox &aabox0, const AABox &aabox1) const;

    bool beyond(double dist, const CoordGroup &group0,
                const CoordGroup &group1) const;

    double minimumDistance(const CoordGroup &group0, const CoordGroup &group1) const;

    double minimumDistance(const AABox &box0, const AABox &box1) const;

    Vector getRandomPoint(const Vector &center, const RanGenerator &generator) const;

	Vector getBoxCenter(const Vector &p) const;
    Vector getBoxCenter(const Vector &p, const Vector &center) const;

    CoordGroup getMinimumImage(const CoordGroup &group, const Vector &center) const;

    CoordGroupArray getMinimumImage(const CoordGroupArray &groups,
                                    const Vector &center,
                                    bool translate_as_one=false) const;

    AABox getMinimumImage(const AABox &aabox, const Vector &center) const;
    
    Vector getMinimumImage(const Vector &point, const Vector &center) const;

    QVector<Vector> getImagesWithin(const Vector &point, const Vector &center, double dist) const;

    QList< boost::tuple<double,CoordGroup> >
               getCopiesWithin(const CoordGroup &group,
                               const CoordGroup &center, double dist) const;

protected:

    Vector wrapDelta(const Vector &v0, const Vector &v1) const;

    CoordGroupArray _pvt_getMinimumImage(
                                const CoordGroupArray &groups,
                                const Vector &point) const;

    static int getWrapVal(double delta, double invlgth, double halflgth);

    /** The lengths of each side of the box */
    Vector boxlength;

    /** Half the box length */
    Vector halflength;

    /** The inverse of the lengths of each side of the box */
    Vector invlength;
};

}

Q_DECLARE_METATYPE(SireVol::PeriodicBox)

SIRE_EXPOSE_CLASS( SireVol::PeriodicBox )

SIRE_END_HEADER

#endif
