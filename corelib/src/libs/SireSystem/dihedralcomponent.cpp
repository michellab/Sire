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

#include "dihedralcomponent.h"
#include "delta.h"

#include "SireVol/space.h"

#include "SireMaths/torsion.h"

#include "SireSystem/system.h"

#include "SireID/index.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireFF;
using namespace SireCAS;
using namespace SireMaths;
using namespace SireVol;
using namespace SireID;
using namespace SireBase;
using namespace SireUnits;
using namespace SireStream;

////////
//////// Implementation of DihedralComponent
////////

static const RegisterMetaType<DihedralComponent> r_dihcomp;

QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds,
                                          const DihedralComponent &dihcomp)
{
    writeHeader(ds, r_dihcomp, 1);
    
    SharedDataStream sds(ds);
    
    sds << dihcomp.p0 << dihcomp.p1 << dihcomp.p2 << dihcomp.p3
        << static_cast<const GeometryComponent&>(dihcomp);
        
    return ds;
}

QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, DihedralComponent &dihcomp)
{
    VersionID v = readHeader(ds, r_dihcomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> dihcomp.p0 >> dihcomp.p1 >> dihcomp.p2 >> dihcomp.p3
            >> static_cast<GeometryComponent&>(dihcomp);
            
        dihcomp.intra_molecule_points = Point::areIntraMoleculePoints(
                                                        dihcomp.p0, dihcomp.p1 ) and
                                        Point::areIntraMoleculePoints(
                                                        dihcomp.p0, dihcomp.p2 ) and
                                        Point::areIntraMoleculePoints(
                                                        dihcomp.p0, dihcomp.p3 );
    }
    else
        throw version_error(v, "1", r_dihcomp, CODELOC);
        
    return ds;
}

/** Null constructor */
DihedralComponent::DihedralComponent() 
               : ConcreteProperty<DihedralComponent,GeometryComponent>(),
                 intra_molecule_points(true)
{}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getPhiSymbol, ("phi") );

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getTheta012Symbol, ("theta012") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getTheta123Symbol, ("theta123") );

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR01Symbol, ("r01") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR12Symbol, ("r12") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR23Symbol, ("r23") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR03Symbol, ("r03") );

/** Return the symbol that represents the dihedral angle */
const Symbol& DihedralComponent::phi()
{
    return *(getPhiSymbol());
}

/** Return the symbol that represents the angle between
    points 0, 1 and 2 */
const Symbol& DihedralComponent::theta012()
{
    return *(getTheta012Symbol());
}

/** Return the symbol that represents the angle between
    points 1, 2 and 3 */
const Symbol& DihedralComponent::theta123()
{
    return *(getTheta123Symbol());
}

/** Return the symbol that represents the 0-1 distance */
const Symbol& DihedralComponent::r01()
{
    return *(getR01Symbol());
}

/** Return the symbol that represents the 1-2 distance */
const Symbol& DihedralComponent::r12()
{
    return *(getR12Symbol());
}

/** Return the symbol that represents the 2-3 distance */
const Symbol& DihedralComponent::r23()
{
    return *(getR23Symbol());
}

/** Return the symbol that represents the 0-3 distance */
const Symbol& DihedralComponent::r03()
{
    return *(getR03Symbol());
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    dihedral between the four points 'point0', 'point1', 'point2' and 'point2' */
DihedralComponent::DihedralComponent(const Symbol &constrained_symbol,
                               const PointRef &point0, const PointRef &point1,
                               const PointRef &point2, const PointRef &point3,
                               const PropertyMap &map)
       : ConcreteProperty<DihedralComponent,GeometryComponent>(constrained_symbol,
                                                               DihedralComponent::phi(),
                                                               map),
         p0(point0), p1(point1), p2(point2), p3(point3)
{
    intra_molecule_points = Point::areIntraMoleculePoints(p0,p1) and
                            Point::areIntraMoleculePoints(p0,p2) and
                            Point::areIntraMoleculePoints(p0,p3);
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    expression based on the dihedral, angles and distances within the four points
    'point0', 'point1', 'point2' and 'point3' */
DihedralComponent::DihedralComponent(const Symbol &constrained_symbol,
                               const PointRef &point0, const PointRef &point1,
                               const PointRef &point2, const PointRef &point3,
                               const Expression &geometry_expression,
                               const PropertyMap &map)
       : ConcreteProperty<DihedralComponent,GeometryComponent>(constrained_symbol,
                                                               geometry_expression,
                                                               map),
         p0(point0), p1(point1), p2(point2), p3(point3)
{
    intra_molecule_points = Point::areIntraMoleculePoints(p0,p1) and
                            Point::areIntraMoleculePoints(p0,p2) and
                            Point::areIntraMoleculePoints(p0,p3);
}
  
/** Copy constructor */
DihedralComponent::DihedralComponent(const DihedralComponent &other)
                  : ConcreteProperty<DihedralComponent,GeometryComponent>(other),
                    p0(other.p0), p1(other.p1), p2(other.p2), p3(other.p3),
                    intra_molecule_points(other.intra_molecule_points)
{}

/** Destructor */
DihedralComponent::~DihedralComponent()
{}

/** Copy assignment operator */
DihedralComponent& DihedralComponent::operator=(const DihedralComponent &other)
{
    if (this != &other)
    {
        GeometryComponent::operator=(other);
        p0 = other.p0;
        p1 = other.p1;
        p2 = other.p2;
        p3 = other.p3;
        intra_molecule_points = other.intra_molecule_points;
    }
    
    return *this;
}

/** Comparison operator */
bool DihedralComponent::operator==(const DihedralComponent &other) const
{
    return this == &other or
           (p0 == other.p0 and p1 == other.p1 and p2 == other.p2 and p3 == other.p3 and
            GeometryComponent::operator==(other));
}

/** Comparison operator */
bool DihedralComponent::operator!=(const DihedralComponent &other) const
{
    return not DihedralComponent::operator==(other);
}

const char* DihedralComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DihedralComponent>() );
}

QString DihedralComponent::toString() const
{
    return QObject::tr("DihedralComponent( %1, %3-%4-%5-%6, %2 )")
                .arg(component().toString(), expression().toString() )
                .arg(p0.read().toString(), p1.read().toString())
                .arg(p2.read().toString(), p3.read().toString());
}

/** Set the space used by the points and distance calculation */
void DihedralComponent::setSpace(const Space &new_space)
{
    if (space().equals(new_space) or intra_molecule_points)
        return;

    DihedralComponent old_state(*this);
    
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
const Point& DihedralComponent::point(int i) const
{
    i = Index(i).map( nPoints() );
    
    switch (i % 4)
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

/** Return the first point between which the dihedral is calculated */
const Point& DihedralComponent::point0() const
{
    return p0.read();
}

/** Return the second point between which the dihedral is calculated */
const Point& DihedralComponent::point1() const
{
    return p1.read();
}

/** Return the third point between which the dihedral is calculated */
const Point& DihedralComponent::point2() const
{
    return p2.read();
}

/** Return the fourth point between which the dihedral is calculated */
const Point& DihedralComponent::point3() const
{
    return p3.read();
}

/** Return the number of points (4) */
int DihedralComponent::nPoints() const
{
    return 4;
}

Torsion DihedralComponent::getTorsion() const
{
    if (intra_molecule_points)
    {
        return Torsion(p0.read().point(), p1.read().point(),
                       p2.read().point(), p3.read().point());
    }
    else
    {
        return Torsion(space().getMinimumImage(p0.read().point(), p2.read().point()),
                       space().getMinimumImage(p1.read().point(), p2.read().point()),
                       p2.read().point(),
                       space().getMinimumImage(p3.read().point(), p2.read().point()));
    }
}

bool DihedralComponent::wouldChange(const Delta &delta, quint32 last_subversion) const
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

Values DihedralComponent::getValues(const System &system)
{
    if (p0.read().wouldUpdate(system))
        p0.edit().update(system);

    if (p1.read().wouldUpdate(system))
        p1.edit().update(system);
        
    if (p2.read().wouldUpdate(system))
        p2.edit().update(system);
        
    if (p3.read().wouldUpdate(system))
        p3.edit().update(system);
        
    Torsion tor = getTorsion();

    Values vals;
    vals.set( phi(), tor.angle().to(radians) );
    vals.set( theta012(), tor.triangle1().angle().to(radians) );
    vals.set( theta123(), tor.triangle2().angle().to(radians) );
    vals.set( r01(), tor.triangle1().line2().length() );
    vals.set( r12(), tor.line12().length() );
    vals.set( r23(), tor.triangle2().line0().length() );
    vals.set( r03(), tor.line03().length() );
        
    return vals;
}
