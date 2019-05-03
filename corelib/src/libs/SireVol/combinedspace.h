/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREVOL_COMBINEDSPACE_H
#define SIREVOL_COMBINEDSPACE_H

#include <QVector>
#include <QList>

#include "space.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class CombinedSpace;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::CombinedSpace&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::CombinedSpace&);

namespace SireVol
{

/** This is a space which is built from the combination of a set
    of sub-spaces. This is useful for systems that are comprised
    of multiple sub-spaces, e.g. Gibbs ensemble simulations,
    or simulations using a combined bound and free binding leg 
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT CombinedSpace 
        : public SireBase::ConcreteProperty<CombinedSpace,Space>
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const CombinedSpace&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, CombinedSpace&);

public:
    typedef QVector<SpacePtr>::const_iterator const_iterator;
    typedef const_iterator iterator;

    CombinedSpace();
    
    CombinedSpace(const Space &space);
    CombinedSpace(const Space &space0, const Space &space1);
    
    CombinedSpace(const QList<SpacePtr> &spaces);
    CombinedSpace(const QVector<SpacePtr> &spaces);
    
    CombinedSpace(const CombinedSpace &other);
    
    ~CombinedSpace();
    
    CombinedSpace& operator=(const CombinedSpace &other);
    
    bool operator==(const CombinedSpace &other) const;
    bool operator!=(const CombinedSpace &other) const;
    
    static const char* typeName();
    
    const Space& operator[](int i) const;
    
    const Space& at(int i) const;
    
    int nSpaces() const;
    int count() const;
    int size();
    
    bool isEmpty() const;

    const_iterator constBegin() const;
    const_iterator begin() const;
    
    const_iterator constEnd() const;
    const_iterator end() const;

    QString toString() const;

    SireUnits::Dimension::Volume volume() const;

    SpacePtr setVolume(SireUnits::Dimension::Volume volume) const;

    double calcDist(const Vector &point0, const Vector &point1) const;
    double calcDist2(const Vector &point0, const Vector &point1) const;

    double calcDist(const CoordGroup &group, DistMatrix &distmat) const;
    double calcDist2(const CoordGroup &group, DistMatrix &distmat) const;

    double calcInvDist(const CoordGroup &group, DistMatrix &distmat) const;
    double calcInvDist2(const CoordGroup &group, DistMatrix &distmat) const;

    double calcDist(const CoordGroup &group1, const CoordGroup &group2,
                    DistMatrix &distmat) const;

    double calcDist2(const CoordGroup &group1, const CoordGroup &group2,
                     DistMatrix &distmat) const;

    double calcDist(const CoordGroup &group, const Vector &point,
                    DistMatrix &distmat) const;
                    
    double calcDist2(const CoordGroup &group, const Vector &point,
                     DistMatrix &distmat) const;

    double calcInvDist(const CoordGroup &group1, const CoordGroup &group2,
                       DistMatrix &distmat) const;

    double calcInvDist2(const CoordGroup &group1, const CoordGroup &group2,
                        DistMatrix &distmat) const;

    DistVector calcDistVector(const Vector &point0, 
                              const Vector &point1) const;

    double calcDistVectors(const CoordGroup &group,
                           DistVectorMatrix &distmat) const;
                                       
    double calcDistVectors(const CoordGroup &group1,
                           const CoordGroup &group2,
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

    bool beyond(double dist, const CoordGroup &group0,
                const CoordGroup &group1) const;

    bool beyond(double dist, const AABox &aabox0, const AABox &aabox1) const;

    double minimumDistance(const AABox &box0, const AABox &box1) const;
    double minimumDistance(const CoordGroup &group0, const CoordGroup &group1) const;
    double minimumDistance(const CoordGroup &group) const;

    bool isPeriodic() const;
    bool isCartesian() const;

    Vector getRandomPoint(const Vector &center, 
                          const RanGenerator &generator) const;

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
               getCopiesWithin(const CoordGroup &group, const CoordGroup &center,
                               double dist) const;

private:
    void assertSingleSpace(const QString &text, const QString &codeloc) const;
    void assertSameSpace(const QString &text, const QString &codeloc) const;

    /** All of the spaces in this combined space */
    QVector<SpacePtr> spces;
};

}

Q_DECLARE_METATYPE( SireVol::CombinedSpace )

SIRE_EXPOSE_CLASS( SireVol::CombinedSpace )

SIRE_END_HEADER

#endif
