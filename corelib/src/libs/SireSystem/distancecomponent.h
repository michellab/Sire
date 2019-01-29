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

#ifndef SIRESYSTEM_DISTANCECOMPONENT_H
#define SIRESYSTEM_DISTANCECOMPONENT_H

#include "geometrycomponent.h"

#include "SireFF/point.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class DistanceComponent;
class DoubleDistanceComponent;
class TripleDistanceComponent;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::DistanceComponent&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::DistanceComponent&);

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::DoubleDistanceComponent&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::DoubleDistanceComponent&);

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::TripleDistanceComponent&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::TripleDistanceComponent&);

namespace SireSystem
{

/** This is a constraint that constrains a symbol to equal the
    value of an expression that involves a distance between atoms
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT DistanceComponent
         : public SireBase::ConcreteProperty<DistanceComponent,GeometryComponent>
{

friend QDataStream& ::operator<<(QDataStream&, const DistanceComponent&);
friend QDataStream& ::operator>>(QDataStream&, DistanceComponent&);

public:
    DistanceComponent();
    DistanceComponent(const SireCAS::Symbol &constrained_symbol,
                      const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                      const PropertyMap &map = PropertyMap());
    DistanceComponent(const SireCAS::Symbol &constrained_symbol,
                      const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                      const SireCAS::Expression &geometry_expression,
                      const PropertyMap &map = PropertyMap());
                      
    DistanceComponent(const DistanceComponent &other);
    
    ~DistanceComponent();
    
    DistanceComponent& operator=(const DistanceComponent &other);
    
    bool operator==(const DistanceComponent &other) const;
    bool operator!=(const DistanceComponent &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const SireFF::Point& point(int i) const;
    
    const SireFF::Point& point0() const;
    const SireFF::Point& point1() const;
    
    int nPoints() const;
    
    static const SireCAS::Symbol& r();
    
protected:
    bool wouldChange(const Delta &delta, quint32 last_subversion) const;
    SireCAS::Values getValues(const System &system);
    
    void setSpace(const SireVol::Space &space);

private:
    double getDistance() const;

    /** The two points between which the distance is calculated */
    SireFF::PointPtr p0, p1;
    
    /** Whether or not the points are within the same molecule */
    bool intra_molecule_points;
};

/** This is a constraint that constrains a symbol to equal the
    value of an expression that involves a distance between 
    two pairs of atoms
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT DoubleDistanceComponent
         : public SireBase::ConcreteProperty<DoubleDistanceComponent,GeometryComponent>
{

friend QDataStream& ::operator<<(QDataStream&, const DoubleDistanceComponent&);
friend QDataStream& ::operator>>(QDataStream&, DoubleDistanceComponent&);

public:
    DoubleDistanceComponent();
    DoubleDistanceComponent(const SireCAS::Symbol &constrained_symbol,
                        const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                        const SireFF::PointRef &point2, const SireFF::PointRef &point3,
                        const PropertyMap &map = PropertyMap());
    DoubleDistanceComponent(const SireCAS::Symbol &constrained_symbol,
                        const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                        const SireFF::PointRef &point2, const SireFF::PointRef &point3,
                        const SireCAS::Expression &geometry_expression,
                        const PropertyMap &map = PropertyMap());
                      
    DoubleDistanceComponent(const DoubleDistanceComponent &other);
    
    ~DoubleDistanceComponent();
    
    DoubleDistanceComponent& operator=(const DoubleDistanceComponent &other);
    
    bool operator==(const DoubleDistanceComponent &other) const;
    bool operator!=(const DoubleDistanceComponent &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const SireFF::Point& point(int i) const;
    
    const SireFF::Point& point0() const;
    const SireFF::Point& point1() const;
    const SireFF::Point& point2() const;
    const SireFF::Point& point3() const;
    
    int nPoints() const;
    
    static const SireCAS::Symbol& r01();
    static const SireCAS::Symbol& r23();
    
protected:
    bool wouldChange(const Delta &delta, quint32 last_subversion) const;
    SireCAS::Values getValues(const System &system);
    
    void setSpace(const SireVol::Space &space);

private:
    double getDistance01() const;
    double getDistance23() const;

    /** The four points between which the distance is calculated */
    SireFF::PointPtr p0, p1, p2, p3;
    
    /** Whether or not the points are within the same molecule */
    bool intra_molecule_points01, intra_molecule_points23;
};

/** This is a constraint that constrains a symbol to equal the
    value of an expression that involves distances between 
    three pairs of atoms
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT TripleDistanceComponent
         : public SireBase::ConcreteProperty<TripleDistanceComponent,GeometryComponent>
{

friend QDataStream& ::operator<<(QDataStream&, const TripleDistanceComponent&);
friend QDataStream& ::operator>>(QDataStream&, TripleDistanceComponent&);

public:
    TripleDistanceComponent();
    TripleDistanceComponent(const SireCAS::Symbol &constrained_symbol,
                       const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                       const SireFF::PointRef &point2, const SireFF::PointRef &point3,
                       const SireFF::PointRef &point4, const SireFF::PointRef &point5,
                       const PropertyMap &map = PropertyMap());
    TripleDistanceComponent(const SireCAS::Symbol &constrained_symbol,
                       const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                       const SireFF::PointRef &point2, const SireFF::PointRef &point3,
                       const SireFF::PointRef &point4, const SireFF::PointRef &point5,
                       const SireCAS::Expression &geometry_expression,
                       const PropertyMap &map = PropertyMap());
                      
    TripleDistanceComponent(const TripleDistanceComponent &other);
    
    ~TripleDistanceComponent();
    
    TripleDistanceComponent& operator=(const TripleDistanceComponent &other);
    
    bool operator==(const TripleDistanceComponent &other) const;
    bool operator!=(const TripleDistanceComponent &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const SireFF::Point& point(int i) const;
    
    const SireFF::Point& point0() const;
    const SireFF::Point& point1() const;
    const SireFF::Point& point2() const;
    const SireFF::Point& point3() const;
    const SireFF::Point& point4() const;
    const SireFF::Point& point5() const;
    
    int nPoints() const;
    
    static const SireCAS::Symbol& r01();
    static const SireCAS::Symbol& r23();
    static const SireCAS::Symbol& r45();
    
protected:
    bool wouldChange(const Delta &delta, quint32 last_subversion) const;
    SireCAS::Values getValues(const System &system);
    
    void setSpace(const SireVol::Space &space);

private:
    double getDistance01() const;
    double getDistance23() const;
    double getDistance45() const;

    /** The points between which the distance is calculated */
    SireFF::PointPtr p0, p1, p2, p3, p4, p5;
    
    /** Whether or not the points are within the same molecule */
    bool intra_molecule_points01, intra_molecule_points23, intra_molecule_points45;
};

}

Q_DECLARE_METATYPE( SireSystem::DistanceComponent )
Q_DECLARE_METATYPE( SireSystem::DoubleDistanceComponent )
Q_DECLARE_METATYPE( SireSystem::TripleDistanceComponent )

SIRE_EXPOSE_CLASS( SireSystem::DistanceComponent )
SIRE_EXPOSE_CLASS( SireSystem::DoubleDistanceComponent )
SIRE_EXPOSE_CLASS( SireSystem::TripleDistanceComponent )

SIRE_END_HEADER

#endif
