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

#ifndef SIRESYSTEM_ANGLECOMPONENT_H
#define SIRESYSTEM_ANGLECOMPONENT_H

#include "geometrycomponent.h"

#include "SireFF/point.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class AngleComponent;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::AngleComponent&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::AngleComponent&);

namespace SireMaths{ class Triangle; }

namespace SireSystem
{

/** This is a constraint that constrains a symbol to equal the
    value of an expression that involves an angle between three points
    or atoms
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT AngleComponent
         : public SireBase::ConcreteProperty<AngleComponent,GeometryComponent>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const AngleComponent&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, AngleComponent&);

public:
    AngleComponent();
    AngleComponent(const SireCAS::Symbol &constrained_symbol,
                   const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                   const SireFF::PointRef &point2,
                   const PropertyMap &map = PropertyMap());
    AngleComponent(const SireCAS::Symbol &constrained_symbol,
                   const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                   const SireFF::PointRef &point2,
                   const SireCAS::Expression &geometry_expression,
                   const PropertyMap &map = PropertyMap());
                      
    AngleComponent(const AngleComponent &other);
    
    ~AngleComponent();
    
    AngleComponent& operator=(const AngleComponent &other);
    
    bool operator==(const AngleComponent &other) const;
    bool operator!=(const AngleComponent &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const SireFF::Point& point(int i) const;
    
    const SireFF::Point& point0() const;
    const SireFF::Point& point1() const;
    const SireFF::Point& point2() const;
    
    int nPoints() const;
    
    static const SireCAS::Symbol& theta();
    static const SireCAS::Symbol& theta102();
    static const SireCAS::Symbol& theta012();
    static const SireCAS::Symbol& theta021();
    
    static const SireCAS::Symbol& r01();
    static const SireCAS::Symbol& r02();
    static const SireCAS::Symbol& r12();
    
protected:
    bool wouldChange(const Delta &delta, quint32 last_subversion) const;
    SireCAS::Values getValues(const System &system);
    
    void setSpace(const SireVol::Space &space);

private:
    SireMaths::Triangle getTriangle() const;

    /** The three points between which the angle is calculated */
    SireFF::PointPtr p0, p1, p2;
    
    /** Whether or not the points are within the same molecule */
    bool intra_molecule_points;
};

}

Q_DECLARE_METATYPE( SireSystem::AngleComponent )

SIRE_EXPOSE_CLASS( SireSystem::AngleComponent )

SIRE_END_HEADER

#endif
