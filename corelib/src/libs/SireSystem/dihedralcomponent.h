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

#ifndef SIRESYSTEM_DIHEDRALCOMPONENT_H
#define SIRESYSTEM_DIHEDRALCOMPONENT_H

#include "geometrycomponent.h"

#include "SireFF/point.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class DihedralComponent;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::DihedralComponent&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::DihedralComponent&);

namespace SireMaths{ class Torsion; }

namespace SireSystem
{

/** This is a constraint that constrains a symbol to equal the
    value of an expression that involves a dihedral between four 
    points or atoms
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT DihedralComponent
         : public SireBase::ConcreteProperty<DihedralComponent,GeometryComponent>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const DihedralComponent&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, DihedralComponent&);

public:
    DihedralComponent();
    DihedralComponent(const SireCAS::Symbol &constrained_symbol,
                      const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                      const SireFF::PointRef &point2, const SireFF::PointRef &point3,
                      const PropertyMap &map = PropertyMap());
    DihedralComponent(const SireCAS::Symbol &constrained_symbol,
                      const SireFF::PointRef &point0, const SireFF::PointRef &point1,
                      const SireFF::PointRef &point2, const SireFF::PointRef &point3,
                      const SireCAS::Expression &geometry_expression,
                      const PropertyMap &map = PropertyMap());
                      
    DihedralComponent(const DihedralComponent &other);
    
    ~DihedralComponent();
    
    DihedralComponent& operator=(const DihedralComponent &other);
    
    bool operator==(const DihedralComponent &other) const;
    bool operator!=(const DihedralComponent &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const SireFF::Point& point(int i) const;
    
    const SireFF::Point& point0() const;
    const SireFF::Point& point1() const;
    const SireFF::Point& point2() const;
    const SireFF::Point& point3() const;
    
    int nPoints() const;
    
    static const SireCAS::Symbol& phi();
    
    static const SireCAS::Symbol& theta012();
    static const SireCAS::Symbol& theta123();
    
    static const SireCAS::Symbol& r01();
    static const SireCAS::Symbol& r12();
    static const SireCAS::Symbol& r23();
    static const SireCAS::Symbol& r03();
    
protected:
    bool wouldChange(const Delta &delta, quint32 last_subversion) const;
    SireCAS::Values getValues(const System &system);
    
    void setSpace(const SireVol::Space &space);

private:
    SireMaths::Torsion getTorsion() const;

    /** The four points between which the dihedral is calculated */
    SireFF::PointPtr p0, p1, p2, p3;
    
    /** Whether or not the points are within the same molecule */
    bool intra_molecule_points;
};

}

Q_DECLARE_METATYPE( SireSystem::DihedralComponent )

SIRE_EXPOSE_CLASS( SireSystem::DihedralComponent )

SIRE_END_HEADER

#endif
