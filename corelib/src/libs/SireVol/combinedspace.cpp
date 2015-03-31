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

#include "combinedspace.h"

#include "coordgroup.h"

#include "SireID/index.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"
#include "SireVol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireVol;
using namespace SireBase;
using namespace SireID;
using namespace SireUnits::Dimension;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<CombinedSpace> r_combinedspace;

/** Serialise to a binary datastream */
QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const CombinedSpace &combined)
{
    writeHeader(ds, r_combinedspace, 1);
    
    SharedDataStream sds(ds);

    sds << combined.spces << static_cast<const Space&>(combined);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, CombinedSpace &combined)
{
    VersionID v = readHeader(ds, r_combinedspace);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> combined.spces >> static_cast<Space&>(combined);
    }
    else
        throw version_error(v, "1", r_combinedspace, CODELOC);

    return ds;
}

/** Construct a default CombinedSpace volume */
CombinedSpace::CombinedSpace() : ConcreteProperty<CombinedSpace,Space>()
{}

/** Construct a combined space from just a single space */    
CombinedSpace::CombinedSpace(const Space &space)
              : ConcreteProperty<CombinedSpace,Space>()
{
    spces.append(space);
    spces.squeeze();
}

/** Construct a combined space from the passed two spaces */
CombinedSpace::CombinedSpace(const Space &space0, const Space &space1)
              : ConcreteProperty<CombinedSpace,Space>()
{
    spces.append(space0);
    spces.append(space1);
    spces.squeeze();
}

/** Construct a combined space from the passed spaces */
CombinedSpace::CombinedSpace(const QList<SpacePtr> &spaces)
              : ConcreteProperty<CombinedSpace,Space>()
{
    if (not spaces.isEmpty())
    {
        spces.reserve( spaces.count() );
        
        for (QList<SpacePtr>::const_iterator it = spaces.constBegin();
             it != spaces.constEnd();
             ++it)
        {
            spces.append(*it);
        }
   
        spces.squeeze();
    }
}

/** Construct a combined space from the passed spaces */
CombinedSpace::CombinedSpace(const QVector<SpacePtr> &spaces)
              : ConcreteProperty<CombinedSpace,Space>(),
                spces(spaces)
{
    spces.squeeze();
}

/** Copy constructor */
CombinedSpace::CombinedSpace(const CombinedSpace &other) 
              : ConcreteProperty<CombinedSpace,Space>(other), spces(other.spces)
{}

/** Destructor */
CombinedSpace::~CombinedSpace()
{}

/** Copy assignment operator */
CombinedSpace& CombinedSpace::operator=(const CombinedSpace &other)
{
    spces = other.spces;
    Space::operator=(other);
    return *this;
}

/** Comparison operator */
bool CombinedSpace::operator==(const CombinedSpace &other) const
{
    return spces == other.spces and Space::operator==(other);
}

/** Comparison operator */
bool CombinedSpace::operator!=(const CombinedSpace &other) const
{
    return spces != other.spces and Space::operator!=(other);
}

/** Return the ith space that makes up this combined space 

    \throw SireError::invalid_index
*/
const Space& CombinedSpace::operator[](int i) const
{
    return spces.at( Index(i).map(spces.count()) );
}

/** Return the ith space that makes up this combined space 

    \throw SireError::invalid_index
*/
const Space& CombinedSpace::at(int i) const
{
    return this->operator[](i);
}

/** Return the number of spaces in this combined space */
int CombinedSpace::nSpaces() const
{
    return spces.count();
}

/** Return the number of spaces in this combined space */
int CombinedSpace::count() const
{
    return this->nSpaces();
}

/** Return the number of spaces in this combined space */
int CombinedSpace::size()
{
    return this->nSpaces();
}

/** Return whether or not this is empty */
bool CombinedSpace::isEmpty() const
{
    return spces.isEmpty();
}

CombinedSpace::const_iterator CombinedSpace::constBegin() const
{
    return spces.constBegin();
}

CombinedSpace::const_iterator CombinedSpace::begin() const
{
    return spces.begin();
}

CombinedSpace::const_iterator CombinedSpace::constEnd() const
{
    return spces.constEnd();
}

CombinedSpace::const_iterator CombinedSpace::end() const
{
    return spces.end();
}

void CombinedSpace::assertSingleSpace(const QString &text, const QString &codeloc) const
{
    if (spces.count() != 1)
        throw SireVol::incompatible_space( QObject::tr(
                "%1 as there is not just a single space. The number of spaces "
                "is equal to %2.").arg(text), codeloc );
}

void CombinedSpace::assertSameSpace(const QString &text, const QString &codeloc) const
{
    if (spces.isEmpty())
        throw SireVol::incompatible_space( QObject::tr(
                "%1 as there are no spaces!").arg(text), codeloc );
    
    else if (spces.count() > 1)
    {
        CombinedSpace::const_iterator it = this->constBegin();
    
        const Space &first_space = it->read();
    
        for (++it; it != this->constEnd(); ++it)
        {
            if ( not first_space.equals(it->read()) )
                throw SireVol::incompatible_space( QObject::tr(
                    "%1 as the spaces are different.")
                        .arg(text), codeloc );
        }
    }
}

/** Return a string representation of this space */
QString CombinedSpace::toString() const
{
    if (spces.isEmpty())
        return QObject::tr( "CombinedSpace{ NULL }" );
    
    else if (spces.count() == 1)
        return spces.at(0).read().toString();
    
    else
    {
        QStringList space_strings;
        
        for (CombinedSpace::const_iterator it = this->constBegin();
             it != this->constEnd();
             ++it)
        {
            space_strings.append( it->read().toString() );
        }
        
        return QString( "CombinedSpace{ %1 }" )
                    .arg(space_strings.join(" : "));
    }
}

/** A CombinedSpace is only periodic of all of the contained
    spaces are periodic */
bool CombinedSpace::isPeriodic() const
{
    if (spces.isEmpty())
        return false;
    
    else if (spces.count() == 1)
        return spces.at(0).read().isPeriodic();
     
    else
    {
        for (CombinedSpace::const_iterator it = this->constBegin();
             it != this->constEnd();
             ++it)
        {
            if (not it->read().isPeriodic())
                return false;
        }
        
        return true;
    }
}

/** A CombinedSpace space is only cartesian if all of the
    sub-spaces are cartesian */
bool CombinedSpace::isCartesian() const
{
    if (spces.isEmpty())
        return false;
    
    else if (spces.count() == 1)
        return spces.at(0).read().isCartesian();
     
    else
    {
        for (CombinedSpace::const_iterator it = this->constBegin();
             it != this->constEnd();
             ++it)
        {
            if (not it->read().isCartesian())
                return false;
        }
        
        return true;
    }
}

/** Returned the combined (summed) volume of all of the spaces */
SireUnits::Dimension::Volume CombinedSpace::volume() const
{
    if (spces.isEmpty())
        return Volume(0);
    
    else if (spces.count() == 1)
        return spces.at(0).read().volume();
     
    else
    {
        Volume total_volume(0);
    
        for (CombinedSpace::const_iterator it = this->constBegin();
             it != this->constEnd();
             ++it)
        {
            total_volume += it->read().volume();
        }
        
        return total_volume;
    }
}

/** Try to set the volume of this combined space - this will only
    work if we are really just one space */
SpacePtr CombinedSpace::setVolume(SireUnits::Dimension::Volume vol) const
{
    this->assertSingleSpace("Cannot set the volume", CODELOC);
    return spces.at(0).read().setVolume(vol);
}

/** Calculate the distance between two points */
double CombinedSpace::calcDist(const Vector &point0, const Vector &point1) const
{
    this->assertSameSpace("Cannot calculate a distance", CODELOC);
    return spces.at(0).read().calcDist(point0, point1);
}

/** Calculate the distance squared between two points */
double CombinedSpace::calcDist2(const Vector &point0, const Vector &point1) const
{
    this->assertSameSpace("Cannot calculate a distance", CODELOC);
    return spces.at(0).read().calcDist2(point0, point1);
}

/** Populate the matrix 'mat' with the distances between all points in
    the group 'group'. Return the shortest distance between points. */
double CombinedSpace::calcDist(const CoordGroup &group, DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDist(group, mat);
}

/** Populate the matrix 'mat' with the distances^2 between all points in
    the group 'group'. Return the shortest distance between points. */
double CombinedSpace::calcDist2(const CoordGroup &group, DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDist2(group, mat);
}

/** Populate the matrix 'mat' with the inverse distances between all points in
    the group 'group'. Return the smallest distance between points. */
double CombinedSpace::calcInvDist(const CoordGroup &group, DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate inverse distances", CODELOC);
    return spces.at(0).read().calcInvDist(group, mat);
}

/** Populate the matrix 'mat' with the inverse distances^2 between all points in
    the group 'group'. Return the smallest distance between points. */
double CombinedSpace::calcInvDist2(const CoordGroup &group, DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate inverse distances", CODELOC);
    return spces.at(0).read().calcInvDist2(group, mat);
}

/** Populate the matrix 'mat' with the distances between all of the
    points of the two CoordGroups. Return the shortest distance between the two
    CoordGroups. */
double CombinedSpace::calcDist(const CoordGroup &group0, const CoordGroup &group1,
                               DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDist(group0, group1, mat);
}

double CombinedSpace::calcDist(const CoordGroup &group, const Vector &point,
                               DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDist(group, point, mat);
}

double CombinedSpace::calcDist2(const CoordGroup &group, const Vector &point,
                               DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDist2(group, point, mat);
}

/** Populate the matrix 'mat' with the distances^2 between all of the
    points of the two CoordGroups. Return the shortest distance between the
    two CoordGroups. */
double CombinedSpace::calcDist2(const CoordGroup &group0, const CoordGroup &group1,
                                DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDist2(group0, group1, mat);
}

/** Populate the matrix 'mat' with the inverse distances between all of the
    points of the two CoordGroups. Return the shortest distance between
    the two CoordGroups. */
double CombinedSpace::calcInvDist(const CoordGroup &group0, const CoordGroup &group1,
                                  DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate inverse distances", CODELOC);
    return spces.at(0).read().calcInvDist(group0, group1, mat);
}

/** Populate the matrix 'mat' with the inverse distances^2 between all of the
    points of the two CoordGroups. Return the shortest distance between
    the two CoordGroups. */
double CombinedSpace::calcInvDist2(const CoordGroup &group0, const CoordGroup &group1,
                                   DistMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate inverse distances", CODELOC);
    return spces.at(0).read().calcInvDist2(group0, group1, mat);
}

/** Calculate the distance vector between two points */
DistVector CombinedSpace::calcDistVector(const Vector &point0, 
                                         const Vector &point1) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDistVector(point0, point1);
}
    
/** Populate the matrix 'distmat' with all of the interpoint distance vectors
    between all points within the CoordGroup. This is *not* a symmetrical matrix,
    as the direction from point A to point B is the negative of the 
    direction from point B to point A. This returns the shortest distance
    between two points in the group (that is not the self-self distance) */
double CombinedSpace::calcDistVectors(const CoordGroup &group, 
                                      DistVectorMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDistVectors(group, mat);
}

/** Populate the matrix 'distmat' between all the points of the two CoordGroups
    'group1' and 'group2' - the returned matrix has the vectors pointing
    from each point in 'group1' to each point in 'group2'. This returns
    the shortest distance between two points in the group */
double CombinedSpace::calcDistVectors(const CoordGroup &group0, const CoordGroup &group1,
                                      DistVectorMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDistVectors(group0, group1, mat);
}

double CombinedSpace::calcDistVectors(const CoordGroup &group, const Vector &point,
                                      DistVectorMatrix &mat) const
{
    this->assertSameSpace("Cannot calculate distances", CODELOC);
    return spces.at(0).read().calcDistVectors(group, point, mat);
}

/** Calculate the angle between the passed three points. This should return
    the acute angle between the points, which should lie between 0 and 180 degrees */
Angle CombinedSpace::calcAngle(const Vector &point0, const Vector &point1,
                               const Vector &point2) const
{
    this->assertSameSpace("Cannot calculate angles", CODELOC);
    return spces.at(0).read().calcAngle(point0, point1, point2);
}

/** Calculate the torsion angle between the passed four points. This should
    return the torsion angle measured clockwise when looking down the 
    torsion from point0-point1-point2-point3. This will lie between 0 and 360 
    degrees */
Angle CombinedSpace::calcDihedral(const Vector &point0, const Vector &point1,
                              const Vector &point2, const Vector &point3) const
{
    this->assertSameSpace("Cannot calculate dihedrals", CODELOC);
    return spces.at(0).read().calcDihedral(point0, point1, point2, point3);
}

/** Return whether or not two groups enclosed by the AABoxes 'aabox0'
    and 'aabox1' are definitely beyond the cutoff distance */
bool CombinedSpace::beyond(double dist, const AABox &aabox0, const AABox &aabox1) const
{
    this->assertSameSpace("Cannot compare positions", CODELOC);
    return spces.at(0).read().beyond(dist, aabox0, aabox1);
}

/** Return whether or not these two groups are definitely beyond the cutoff distance. */
bool CombinedSpace::beyond(double dist, const CoordGroup &group0,
                           const CoordGroup &group1) const
{
    this->assertSameSpace("Cannot compare positions", CODELOC);
    return spces.at(0).read().beyond(dist, group0, group1);
}

/** Return the minimum distance between the points in 'group0' and 'group1'. */
double CombinedSpace::minimumDistance(const CoordGroup &group0,
                                      const CoordGroup &group1) const
{
    this->assertSameSpace("Cannot calculate minimum distances", CODELOC);
    return spces.at(0).read().minimumDistance(group0, group1);
}

/** Return the minimum distance between the boxes 'box0' and 'box1'. */
double CombinedSpace::minimumDistance(const AABox &box0, const AABox &box1) const
{
    this->assertSameSpace("Cannot calculate minimum distances", CODELOC);
    return spces.at(0).read().minimumDistance(box0, box1);
}

/** Return the minimum distance between points within the group 'group'. */
double CombinedSpace::minimumDistance(const CoordGroup &group) const
{
    this->assertSameSpace("Cannot calculate minimum distances", CODELOC);
    return spces.at(0).read().minimumDistance(group);
}

/** Return the minimum image copy of 'group' with respect to 'center'.
    In this case, as this is not a periodic space, this just returns
    'group' */
CoordGroup CombinedSpace::getMinimumImage(const CoordGroup &group, 
                                          const Vector &center) const
{
    this->assertSameSpace("Cannot get minimum images", CODELOC);
    return spces.at(0).read().getMinimumImage(group, center);
}

/** Return the minimum image copy of 'groups' with respect to 'center'.
    In this case, as this is not a periodic space, this just returns
    'groups' */
CoordGroupArray CombinedSpace::getMinimumImage(const CoordGroupArray &groups,
                                               const Vector &center, bool asone) const
{
    this->assertSameSpace("Cannot get minimum images", CODELOC);
    return spces.at(0).read().getMinimumImage(groups, center, asone);
}

/** A cartesian space is not periodic, so this just returns the input 'aabox' */
AABox CombinedSpace::getMinimumImage(const AABox &aabox, const Vector &center) const
{
    this->assertSameSpace("Cannot get minimum images", CODELOC);
    return spces.at(0).read().getMinimumImage(aabox, center);
}

/** A cartesian space is not periodic, so this just returns the input 'point' */
Vector CombinedSpace::getMinimumImage(const Vector &point, const Vector &center) const
{
    this->assertSameSpace("Cannot get minimum images", CODELOC);
    return spces.at(0).read().getMinimumImage(point, center);
}

QVector<Vector> CombinedSpace::getImagesWithin(const Vector &point, const Vector &center,
                                               double dist) const
{
    this->assertSameSpace("Cannot get minimum images", CODELOC);
    return spces.at(0).read().getImagesWithin(point, center, dist);
}

/** Return a random point within the spaces used in this combined space */
Vector CombinedSpace::getRandomPoint(const Vector &center,
                                     const RanGenerator &generator) const
{
    if (spces.isEmpty())
        return center;
    
    else if (spces.count() == 1)
    {
        return spces.at(0).read().getRandomPoint(center, generator);
    }
    else
    {
        //choose a space at random
        int idx = generator.randInt( spces.count() - 1 );
        
        return spces.at(idx).read().getRandomPoint(center, generator);
    }
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at the origin */
Vector CombinedSpace::getBoxCenter(const Vector &p) const
{
	this->assertSameSpace("Cannot get the box center", CODELOC);
    return spces.at(0).read().getBoxCenter(p);
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at 'center' */
Vector CombinedSpace::getBoxCenter(const Vector &p, const Vector &center) const
{
	this->assertSameSpace("Cannot get the box center", CODELOC);
    return spces.at(0).read().getBoxCenter(p, center);
}

/** Return a list of copies of CoordGroup 'group' that are within
    'distance' of the CoordGroup 'center', translating 'group' so that
    it has the right coordinates to be around 'center'. As this is not
    a periodic space, this will merely return a copy of 'group' if
    it is within the specified distance. */
QList< tuple<double,CoordGroup> >
CombinedSpace::getCopiesWithin(const CoordGroup &group, const CoordGroup &center,
                               double dist) const
{
    this->assertSameSpace("Cannot get close copies", CODELOC);
    return spces.at(0).read().getCopiesWithin(group, center, dist);
}

const char* CombinedSpace::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CombinedSpace>() );
}
