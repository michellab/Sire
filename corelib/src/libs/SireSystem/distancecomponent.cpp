/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "distancecomponent.h"
#include "delta.h"

#include "SireVol/space.h"

#include "SireSystem/system.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireFF;
using namespace SireCAS;
using namespace SireVol;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

////////
//////// Implementation of DistanceComponent
////////

static const RegisterMetaType<DistanceComponent> r_distcomp;

QDataStream &operator<<(QDataStream &ds,
                                          const DistanceComponent &distcomp)
{
    writeHeader(ds, r_distcomp, 1);
    
    SharedDataStream sds(ds);
    
    sds << distcomp.p0 << distcomp.p1
        << static_cast<const GeometryComponent&>(distcomp);
        
    return ds;
}

QDataStream &operator>>(QDataStream &ds, DistanceComponent &distcomp)
{
    VersionID v = readHeader(ds, r_distcomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> distcomp.p0 >> distcomp.p1
            >> static_cast<GeometryComponent&>(distcomp);
            
        distcomp.intra_molecule_points = Point::areIntraMoleculePoints(
                                                        distcomp.p0, distcomp.p1 );
    }
    else
        throw version_error(v, "1", r_distcomp, CODELOC);
        
    return ds;
}

/** Null constructor */
DistanceComponent::DistanceComponent() 
                  : ConcreteProperty<DistanceComponent,GeometryComponent>(),
                    intra_molecule_points(true)
{}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getRSymbol, ("r") );

/** Return the symbol that represents the distance between the
    two points ("r") */
const Symbol& DistanceComponent::r()
{
    return *(getRSymbol());
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    distance between the two points 'point0' and 'point1' */
DistanceComponent::DistanceComponent(const Symbol &constrained_symbol,
                                     const PointRef &point0, const PointRef &point1,
                                     const PropertyMap &map)
       : ConcreteProperty<DistanceComponent,GeometryComponent>(constrained_symbol,
                                                               DistanceComponent::r(),
                                                               map),
         p0(point0), p1(point1)
{
    intra_molecule_points = Point::areIntraMoleculePoints(p0,p1);
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    expression based on the distance between the two points
    'point0' and 'point1' */
DistanceComponent::DistanceComponent(const Symbol &constrained_symbol,
                                     const PointRef &point0, const PointRef &point1,
                                     const Expression &geometry_expression,
                                     const PropertyMap &map)
       : ConcreteProperty<DistanceComponent,GeometryComponent>(constrained_symbol,
                                                               geometry_expression,
                                                               map),
         p0(point0), p1(point1)
{
    intra_molecule_points = Point::areIntraMoleculePoints(p0,p1);
}
  
/** Copy constructor */
DistanceComponent::DistanceComponent(const DistanceComponent &other)
                  : ConcreteProperty<DistanceComponent,GeometryComponent>(other),
                    p0(other.p0), p1(other.p1),
                    intra_molecule_points(other.intra_molecule_points)
{}

/** Destructor */
DistanceComponent::~DistanceComponent()
{}

/** Copy assignment operator */
DistanceComponent& DistanceComponent::operator=(const DistanceComponent &other)
{
    if (this != &other)
    {
        GeometryComponent::operator=(other);
        p0 = other.p0;
        p1 = other.p1;
        intra_molecule_points = other.intra_molecule_points;
    }
    
    return *this;
}

/** Comparison operator */
bool DistanceComponent::operator==(const DistanceComponent &other) const
{
    return this == &other or
           (p0 == other.p0 and p1 == other.p1 and 
            GeometryComponent::operator==(other));
}

/** Comparison operator */
bool DistanceComponent::operator!=(const DistanceComponent &other) const
{
    return not DistanceComponent::operator==(other);
}

const char* DistanceComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DistanceComponent>() );
}

QString DistanceComponent::toString() const
{
    return QObject::tr("DistanceComponent( %1, %3-%4, %2 )")
                .arg(component().toString(), expression().toString() )
                .arg(p0.read().toString(), p1.read().toString());
}

/** Set the space used by the points and distance calculation */
void DistanceComponent::setSpace(const Space &new_space)
{
    if (space().equals(new_space) or intra_molecule_points)
        return;

    DistanceComponent old_state(*this);
    
    try
    {
        p0.edit().setSpace(new_space);
        p1.edit().setSpace(new_space);
        GeometryComponent::setSpace(new_space);
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Return the ith point

    \throw SireError::invalid_index
*/
const Point& DistanceComponent::point(int i) const
{
    i = Index(i).map( nPoints() );
    
    if (i == 0)
        return p0.read();
    else
        return p1.read();
}

/** Return the first point between which the distance is calculated */
const Point& DistanceComponent::point0() const
{
    return p0.read();
}

/** Return the second point between which the distance is calculated */
const Point& DistanceComponent::point1() const
{
    return p1.read();
}

/** Return the number of points (2) */
int DistanceComponent::nPoints() const
{
    return 2;
}

double DistanceComponent::getDistance() const
{
    if (intra_molecule_points)
        return Vector::distance(p0.read().point(), p1.read().point());
    else
        return space().calcDist(p0.read().point(), p1.read().point());
}

bool DistanceComponent::wouldChange(const Delta &delta, quint32 last_subversion) const
{
    if (delta.hasMoleculeChangeSince(last_subversion))
    {
        return delta.sinceChanged(p0.read(), last_subversion) or
               delta.sinceChanged(p1.read(), last_subversion);
    }
    else
        return false;
}

Values DistanceComponent::getValues(const System &system)
{
    if (p0.read().wouldUpdate(system))
        p0.edit().update(system);

    if (p1.read().wouldUpdate(system))
        p1.edit().update(system);
        
    return r() == getDistance();
}

////////
//////// Implementation of DoubleDistanceComponent
////////

static const RegisterMetaType<DoubleDistanceComponent> r_dist2comp;

QDataStream &operator<<(QDataStream &ds,
                                          const DoubleDistanceComponent &dist2comp)
{
    writeHeader(ds, r_dist2comp, 1);
    
    SharedDataStream sds(ds);
    
    sds << dist2comp.p0 << dist2comp.p1 << dist2comp.p2 << dist2comp.p3
        << static_cast<const GeometryComponent&>(dist2comp);
        
    return ds;
}

QDataStream &operator>>(QDataStream &ds, 
                                          DoubleDistanceComponent &dist2comp)
{
    VersionID v = readHeader(ds, r_dist2comp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> dist2comp.p0 >> dist2comp.p1 >> dist2comp.p2 >> dist2comp.p3
            >> static_cast<GeometryComponent&>(dist2comp);
            
        dist2comp.intra_molecule_points01 = Point::areIntraMoleculePoints(
                                                        dist2comp.p0, dist2comp.p1 );
        
        dist2comp.intra_molecule_points23 = Point::areIntraMoleculePoints(
                                                        dist2comp.p2, dist2comp.p3 );
        
    }
    else
        throw version_error(v, "1", r_dist2comp, CODELOC);
        
    return ds;
}

/** Null constructor */
DoubleDistanceComponent::DoubleDistanceComponent() 
                        : ConcreteProperty<DoubleDistanceComponent,GeometryComponent>(),
                          intra_molecule_points01(true), intra_molecule_points23(true)
{}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR01Symbol, ("r01") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR23Symbol, ("r23") );

/** Return the symbol that represents the distance between the
    points .point0() and .point1() ("r01") */
const Symbol& DoubleDistanceComponent::r01()
{
    return *(getR01Symbol());
}

/** Return the symbol that represents the distance between the
    points .point2() and .point2() ("r23") */
const Symbol& DoubleDistanceComponent::r23()
{
    return *(getR23Symbol());
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    distance between the two points 'point0' and 'point1' */
DoubleDistanceComponent::DoubleDistanceComponent(const Symbol &constrained_symbol,
                                      const PointRef &point0, const PointRef &point1,
                                      const PointRef &point2, const PointRef &point3,
                                      const PropertyMap &map)
       : ConcreteProperty<DoubleDistanceComponent,GeometryComponent>(constrained_symbol,
                        DoubleDistanceComponent::r01() + DoubleDistanceComponent::r23(),
                        map),
         p0(point0), p1(point1), p2(point2), p3(point3)
{
    intra_molecule_points01 = Point::areIntraMoleculePoints(p0,p1);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p2,p3);
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    expression based on the distance between the two points
    'point0' and 'point1' */
DoubleDistanceComponent::DoubleDistanceComponent(const Symbol &constrained_symbol,
                                     const PointRef &point0, const PointRef &point1,
                                     const PointRef &point2, const PointRef &point3,
                                     const Expression &geometry_expression,
                                     const PropertyMap &map)
       : ConcreteProperty<DoubleDistanceComponent,GeometryComponent>(constrained_symbol,
                                                               geometry_expression, map),
         p0(point0), p1(point1), p2(point2), p3(point3)
{
    intra_molecule_points01 = Point::areIntraMoleculePoints(p0,p1);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p2,p3);
}
  
/** Copy constructor */
DoubleDistanceComponent::DoubleDistanceComponent(const DoubleDistanceComponent &other)
                  : ConcreteProperty<DoubleDistanceComponent,GeometryComponent>(other),
                    p0(other.p0), p1(other.p1), p2(other.p2), p3(other.p3),
                    intra_molecule_points01(other.intra_molecule_points01),
                    intra_molecule_points23(other.intra_molecule_points23)
{}

/** Destructor */
DoubleDistanceComponent::~DoubleDistanceComponent()
{}

/** Copy assignment operator */
DoubleDistanceComponent& 
DoubleDistanceComponent::operator=(const DoubleDistanceComponent &other)
{
    if (this != &other)
    {
        GeometryComponent::operator=(other);
        p0 = other.p0;
        p1 = other.p1;
        p2 = other.p2;
        p3 = other.p3;
        intra_molecule_points01 = other.intra_molecule_points01;
        intra_molecule_points23 = other.intra_molecule_points23;
    }
    
    return *this;
}

/** Comparison operator */
bool DoubleDistanceComponent::operator==(const DoubleDistanceComponent &other) const
{
    return this == &other or
           (p0 == other.p0 and p1 == other.p1 and 
            p2 == other.p2 and p3 == other.p3 and
            GeometryComponent::operator==(other));
}

/** Comparison operator */
bool DoubleDistanceComponent::operator!=(const DoubleDistanceComponent &other) const
{
    return not DoubleDistanceComponent::operator==(other);
}

const char* DoubleDistanceComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DoubleDistanceComponent>() );
}

QString DoubleDistanceComponent::toString() const
{
    return QObject::tr("DoubleDistanceComponent( %1, %3-%4 and %5-%6, %2 )")
                .arg(component().toString(), expression().toString() )
                .arg(p0.read().toString(), p1.read().toString())
                .arg(p2.read().toString(), p3.read().toString());
}

/** Set the space used by the points and distance calculation */
void DoubleDistanceComponent::setSpace(const Space &new_space)
{
    if (space().equals(new_space) or 
        (intra_molecule_points01 and intra_molecule_points23) )
    {
        return;
    }

    DoubleDistanceComponent old_state(*this);
    
    try
    {
        p0.edit().setSpace(new_space);
        p1.edit().setSpace(new_space);
        p2.edit().setSpace(new_space);
        p3.edit().setSpace(new_space);
        GeometryComponent::setSpace(new_space);
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Return the ith point

    \throw SireError::invalid_index
*/
const Point& DoubleDistanceComponent::point(int i) const
{
    i = Index(i).map( nPoints() );
    
    switch (i)
    {
    case 0:
        return p0.read();
    case 1:
        return p1.read();
    case 2:
        return p2.read();
    default:
        return p3.read();
    }
}

/** Return the first point between which the first distance is calculated */
const Point& DoubleDistanceComponent::point0() const
{
    return p0.read();
}

/** Return the second point between which the first distance is calculated */
const Point& DoubleDistanceComponent::point1() const
{
    return p1.read();
}

/** Return the first point between which the second distance is calculated */
const Point& DoubleDistanceComponent::point2() const
{
    return p2.read();
}

/** Return the second point between which the second distance is calculated */
const Point& DoubleDistanceComponent::point3() const
{
    return p3.read();
}

/** Return the number of points (4) */
int DoubleDistanceComponent::nPoints() const
{
    return 4;
}

double DoubleDistanceComponent::getDistance01() const
{
    if (intra_molecule_points01)
        return Vector::distance(p0.read().point(), p1.read().point());
    else
        return space().calcDist(p0.read().point(), p1.read().point());
}

double DoubleDistanceComponent::getDistance23() const
{
    if (intra_molecule_points23)
        return Vector::distance(p2.read().point(), p3.read().point());
    else
        return space().calcDist(p2.read().point(), p3.read().point());
}

bool DoubleDistanceComponent::wouldChange(const Delta &delta, 
                                          quint32 last_subversion) const
{
    if (delta.hasMoleculeChangeSince(last_subversion))
    {
        return delta.sinceChanged(p0.read(), last_subversion) or
               delta.sinceChanged(p1.read(), last_subversion) or
               delta.sinceChanged(p2.read(), last_subversion) or
               delta.sinceChanged(p3.read(), last_subversion);
    }
    else
        return false;
}

Values DoubleDistanceComponent::getValues(const System &system)
{
    if (p0.read().wouldUpdate(system))
        p0.edit().update(system);

    if (p1.read().wouldUpdate(system))
        p1.edit().update(system);

    if (p2.read().wouldUpdate(system))
        p2.edit().update(system);

    if (p3.read().wouldUpdate(system))
        p3.edit().update(system);
        
    Values vals;
    
    vals.set( r01(), getDistance01() );
    vals.set( r23(), getDistance23() );
        
    return vals;
}

////////
//////// Implementation of TripleDistanceComponent
////////

static const RegisterMetaType<TripleDistanceComponent> r_dist3comp;

QDataStream &operator<<(QDataStream &ds,
                                          const TripleDistanceComponent &dist3comp)
{
    writeHeader(ds, r_dist3comp, 1);
    
    SharedDataStream sds(ds);
    
    sds << dist3comp.p0 << dist3comp.p1 << dist3comp.p2 << dist3comp.p3
        << dist3comp.p4 << dist3comp.p5
        << static_cast<const GeometryComponent&>(dist3comp);
        
    return ds;
}

QDataStream &operator>>(QDataStream &ds, 
                                          TripleDistanceComponent &dist3comp)
{
    VersionID v = readHeader(ds, r_dist3comp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> dist3comp.p0 >> dist3comp.p1 >> dist3comp.p2 >> dist3comp.p3
            >> dist3comp.p4 >> dist3comp.p5
            >> static_cast<GeometryComponent&>(dist3comp);
            
        dist3comp.intra_molecule_points01 = Point::areIntraMoleculePoints(
                                                        dist3comp.p0, dist3comp.p1 );
        
        dist3comp.intra_molecule_points23 = Point::areIntraMoleculePoints(
                                                        dist3comp.p2, dist3comp.p3 );
        
        dist3comp.intra_molecule_points45 = Point::areIntraMoleculePoints(
                                                        dist3comp.p4, dist3comp.p5 );
    }
    else
        throw version_error(v, "1", r_dist3comp, CODELOC);
        
    return ds;
}

/** Null constructor */
TripleDistanceComponent::TripleDistanceComponent() 
                        : ConcreteProperty<TripleDistanceComponent,GeometryComponent>(),
                          intra_molecule_points01(true), intra_molecule_points23(true),
                          intra_molecule_points45(true)
{}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR45Symbol, ("r45") );

/** Return the symbol that represents the distance between the
    points .point0() and .point1() ("r01") */
const Symbol& TripleDistanceComponent::r01()
{
    return *(getR01Symbol());
}

/** Return the symbol that represents the distance between the
    points .point2() and .point3() ("r23") */
const Symbol& TripleDistanceComponent::r23()
{
    return *(getR23Symbol());
}

/** Return the symbol that represents the distance between the
    points .point4() and .point5() ("r45") */
const Symbol& TripleDistanceComponent::r45()
{
    return *(getR45Symbol());
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    distance between the two points 'point0' and 'point1' */
TripleDistanceComponent::TripleDistanceComponent(const Symbol &constrained_symbol,
                                      const PointRef &point0, const PointRef &point1,
                                      const PointRef &point2, const PointRef &point3,
                                      const PointRef &point4, const PointRef &point5,
                                      const PropertyMap &map)
       : ConcreteProperty<TripleDistanceComponent,GeometryComponent>(constrained_symbol,
                        TripleDistanceComponent::r01() + TripleDistanceComponent::r23(),
                        map),
         p0(point0), p1(point1), p2(point2), p3(point3),
         p4(point4), p5(point5)
{
    intra_molecule_points01 = Point::areIntraMoleculePoints(p0,p1);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p2,p3);
    intra_molecule_points45 = Point::areIntraMoleculePoints(p4,p5);
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    expression based on the distance between the two points
    'point0' and 'point1' */
TripleDistanceComponent::TripleDistanceComponent(const Symbol &constrained_symbol,
                                     const PointRef &point0, const PointRef &point1,
                                     const PointRef &point2, const PointRef &point3,
                                     const PointRef &point4, const PointRef &point5,
                                     const Expression &geometry_expression,
                                     const PropertyMap &map)
       : ConcreteProperty<TripleDistanceComponent,GeometryComponent>(constrained_symbol,
                                                               geometry_expression, map),
         p0(point0), p1(point1), p2(point2), p3(point3),
         p4(point4), p5(point5)
{
    intra_molecule_points01 = Point::areIntraMoleculePoints(p0,p1);
    intra_molecule_points23 = Point::areIntraMoleculePoints(p2,p3);
    intra_molecule_points45 = Point::areIntraMoleculePoints(p4,p5);
}
  
/** Copy constructor */
TripleDistanceComponent::TripleDistanceComponent(const TripleDistanceComponent &other)
                  : ConcreteProperty<TripleDistanceComponent,GeometryComponent>(other),
                    p0(other.p0), p1(other.p1), p2(other.p2), p3(other.p3),
                    p4(other.p4), p5(other.p5),
                    intra_molecule_points01(other.intra_molecule_points01),
                    intra_molecule_points23(other.intra_molecule_points23),
                    intra_molecule_points45(other.intra_molecule_points45)
{}

/** Destructor */
TripleDistanceComponent::~TripleDistanceComponent()
{}

/** Copy assignment operator */
TripleDistanceComponent& 
TripleDistanceComponent::operator=(const TripleDistanceComponent &other)
{
    if (this != &other)
    {
        GeometryComponent::operator=(other);
        p0 = other.p0;
        p1 = other.p1;
        p2 = other.p2;
        p3 = other.p3;
        p4 = other.p4;
        p5 = other.p5;
        intra_molecule_points01 = other.intra_molecule_points01;
        intra_molecule_points23 = other.intra_molecule_points23;
        intra_molecule_points45 = other.intra_molecule_points45;
    }
    
    return *this;
}

/** Comparison operator */
bool TripleDistanceComponent::operator==(const TripleDistanceComponent &other) const
{
    return this == &other or
           (p0 == other.p0 and p1 == other.p1 and 
            p2 == other.p2 and p3 == other.p3 and
            p4 == other.p4 and p5 == other.p5 and
            GeometryComponent::operator==(other));
}

/** Comparison operator */
bool TripleDistanceComponent::operator!=(const TripleDistanceComponent &other) const
{
    return not TripleDistanceComponent::operator==(other);
}

const char* TripleDistanceComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TripleDistanceComponent>() );
}

QString TripleDistanceComponent::toString() const
{
    return QObject::tr("TripleDistanceComponent( %1, %3-%4, %5-%6 and %7-%8, %2 )")
                .arg(component().toString(), expression().toString() )
                .arg(p0.read().toString(), p1.read().toString())
                .arg(p2.read().toString(), p3.read().toString())
                .arg(p4.read().toString(), p5.read().toString());
}

/** Set the space used by the points and distance calculation */
void TripleDistanceComponent::setSpace(const Space &new_space)
{
    if (space().equals(new_space) or 
        (intra_molecule_points01 and intra_molecule_points23 and intra_molecule_points45))
    {
        return;
    }

    TripleDistanceComponent old_state(*this);
    
    try
    {
        p0.edit().setSpace(new_space);
        p1.edit().setSpace(new_space);
        p2.edit().setSpace(new_space);
        p3.edit().setSpace(new_space);
        p4.edit().setSpace(new_space);
        p5.edit().setSpace(new_space);
        GeometryComponent::setSpace(new_space);
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Return the ith point

    \throw SireError::invalid_index
*/
const Point& TripleDistanceComponent::point(int i) const
{
    i = Index(i).map( nPoints() );
    
    switch (i)
    {
    case 0:
        return p0.read();
    case 1:
        return p1.read();
    case 2:
        return p2.read();
    case 3:
        return p3.read();
    case 4:
        return p4.read();
    default:
        return p5.read();
    }
}

/** Return the first point between which the first distance is calculated */
const Point& TripleDistanceComponent::point0() const
{
    return p0.read();
}

/** Return the second point between which the first distance is calculated */
const Point& TripleDistanceComponent::point1() const
{
    return p1.read();
}

/** Return the first point between which the second distance is calculated */
const Point& TripleDistanceComponent::point2() const
{
    return p2.read();
}

/** Return the second point between which the second distance is calculated */
const Point& TripleDistanceComponent::point3() const
{
    return p3.read();
}

/** Return the first point between which the third distance is calculated */
const Point& TripleDistanceComponent::point4() const
{
    return p4.read();
}

/** Return the second point between which the third distance is calculated */
const Point& TripleDistanceComponent::point5() const
{
    return p5.read();
}

/** Return the number of points (6) */
int TripleDistanceComponent::nPoints() const
{
    return 6;
}

double TripleDistanceComponent::getDistance01() const
{
    if (intra_molecule_points01)
        return Vector::distance(p0.read().point(), p1.read().point());
    else
        return space().calcDist(p0.read().point(), p1.read().point());
}

double TripleDistanceComponent::getDistance23() const
{
    if (intra_molecule_points23)
        return Vector::distance(p2.read().point(), p3.read().point());
    else
        return space().calcDist(p2.read().point(), p3.read().point());
}

double TripleDistanceComponent::getDistance45() const
{
    if (intra_molecule_points45)
        return Vector::distance(p4.read().point(), p5.read().point());
    else
        return space().calcDist(p4.read().point(), p5.read().point());
}

bool TripleDistanceComponent::wouldChange(const Delta &delta,
                                          quint32 last_subversion) const
{
    if (delta.hasMoleculeChangeSince(last_subversion))
    {
        return delta.sinceChanged(p0.read(), last_subversion) or
               delta.sinceChanged(p1.read(), last_subversion) or
               delta.sinceChanged(p2.read(), last_subversion) or
               delta.sinceChanged(p3.read(), last_subversion);
    }
    else
        return false;
}

Values TripleDistanceComponent::getValues(const System &system)
{
    if (p0.read().wouldUpdate(system))
        p0.edit().update(system);

    if (p1.read().wouldUpdate(system))
        p1.edit().update(system);

    if (p2.read().wouldUpdate(system))
        p2.edit().update(system);

    if (p3.read().wouldUpdate(system))
        p3.edit().update(system);

    if (p4.read().wouldUpdate(system))
        p4.edit().update(system);

    if (p5.read().wouldUpdate(system))
        p5.edit().update(system);
        
    Values vals;
    
    vals.set( r01(), getDistance01() );
    vals.set( r23(), getDistance23() );
    vals.set( r45(), getDistance45() );
        
    return vals;
}
