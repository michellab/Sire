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

#ifndef SIREVOL_CARTESIAN_H
#define SIREVOL_CARTESIAN_H

#include "space.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class Cartesian;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::Cartesian&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::Cartesian&);

namespace SireVol
{

/**
This class overloads SimVolume to provide an infinite Cartesian
(3-dimensional, orthoganol dimensions) volume. This corresponds to
a traditional gas-phase or no-boundary system.

@author Christopher Woods
*/
class SIREVOL_EXPORT Cartesian
          : public SireBase::ConcreteProperty<Cartesian,Space>
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Cartesian&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, Cartesian&);

public:
    Cartesian();
    Cartesian(const Cartesian &other);

    virtual ~Cartesian();

    Cartesian& operator=(const Cartesian &other);

    bool operator==(const Cartesian &other) const;
    bool operator!=(const Cartesian &other) const;

    static const char* typeName();

    QString toString() const;

    bool isPeriodic() const;
    bool isCartesian() const;

    SireUnits::Dimension::Volume volume() const;
    SpacePtr setVolume(SireUnits::Dimension::Volume volume) const;

    double calcDist(const Vector &point0, const Vector &point1) const;
    double calcDist2(const Vector &point0, const Vector &point1) const;

    double calcDist(const CoordGroup &group, DistMatrix &mat) const;
    double calcDist2(const CoordGroup &group, DistMatrix &mat) const;
    double calcInvDist(const CoordGroup &group, DistMatrix &mat) const;
    double calcInvDist2(const CoordGroup &group, DistMatrix &mat) const;

    double calcDist(const CoordGroup &group1, const CoordGroup &group2,
                    DistMatrix &mat) const;

    double calcDist(const CoordGroup &group, const Vector &point,
                    DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group, const Vector &point,
                     DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group1, const CoordGroup &group2,
                     DistMatrix &mat) const;

    double calcInvDist(const CoordGroup &group1, const CoordGroup &group2,
                       DistMatrix &mat) const;

    double calcInvDist2(const CoordGroup &group1, const CoordGroup &group2,
                        DistMatrix &mat) const;

    DistVector calcDistVector(const Vector &point0, const Vector &point1) const;

    double calcDistVectors(const CoordGroup &group, DistVectorMatrix &distmat) const;
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

    double minimumDistance(const CoordGroup &group) const;

    double minimumDistance(const AABox &box0, const AABox &box1) const;

    double minimumDistance(const Vector &p, const AABox &box) const;

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
};

}

Q_DECLARE_METATYPE(SireVol::Cartesian)

SIRE_EXPOSE_CLASS( SireVol::Cartesian )

SIRE_END_HEADER

#endif
