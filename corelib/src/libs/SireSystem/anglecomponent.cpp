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

#include "anglecomponent.h"
#include "delta.h"

#include "SireVol/space.h"

#include "SireMaths/triangle.h"

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
//////// Implementation of AngleComponent
////////

static const RegisterMetaType<AngleComponent> r_angcomp;

QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds,
                                          const AngleComponent &angcomp)
{
    writeHeader(ds, r_angcomp, 1);
    
    SharedDataStream sds(ds);
    
    sds << angcomp.p0 << angcomp.p1 << angcomp.p2
        << static_cast<const GeometryComponent&>(angcomp);
        
    return ds;
}

QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, AngleComponent &angcomp)
{
    VersionID v = readHeader(ds, r_angcomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> angcomp.p0 >> angcomp.p1 >> angcomp.p2
            >> static_cast<GeometryComponent&>(angcomp);
            
        angcomp.intra_molecule_points = Point::areIntraMoleculePoints(
                                                        angcomp.p0, angcomp.p1 ) and
                                        Point::areIntraMoleculePoints(
                                                        angcomp.p0, angcomp.p2 );
    }
    else
        throw version_error(v, "1", r_angcomp, CODELOC);
        
    return ds;
}

/** Null constructor */
AngleComponent::AngleComponent() 
               : ConcreteProperty<AngleComponent,GeometryComponent>(),
                 intra_molecule_points(true)
{}

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getThetaSymbol, ("theta") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getTheta012Symbol, ("theta012") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getTheta102Symbol, ("theta102") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getTheta021Symbol, ("theta021") );

Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR01Symbol, ("r01") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR02Symbol, ("r02") );
Q_GLOBAL_STATIC_WITH_ARGS( Symbol, getR12Symbol, ("r12") );

/** Return the symbol that represents the central angle between
    the three points*/
const Symbol& AngleComponent::theta()
{
    return *(getThetaSymbol());
}

/** Return the symbol that represents the angle 012 between
    the three points*/
const Symbol& AngleComponent::theta012()
{
    return *(getTheta012Symbol());
}

/** Return the symbol that represents the angle 102 between
    the three points*/
const Symbol& AngleComponent::theta102()
{
    return *(getTheta102Symbol());
}

/** Return the symbol that represents the angle 021 between
    the three points*/
const Symbol& AngleComponent::theta021()
{
    return *(getTheta021Symbol());
}

/** Return the symbol that represents the 0-1 distance */
const Symbol& AngleComponent::r01()
{
    return *(getR01Symbol());
}

/** Return the symbol that represents the 0-2 distance */
const Symbol& AngleComponent::r02()
{
    return *(getR02Symbol());
}

/** Return the symbol that represents the 1-2 distance */
const Symbol& AngleComponent::r12()
{
    return *(getR12Symbol());
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    angle between the three points 'point0', 'point1' and 'point2' */
AngleComponent::AngleComponent(const Symbol &constrained_symbol,
                               const PointRef &point0, const PointRef &point1,
                               const PointRef &point2,
                               const PropertyMap &map)
       : ConcreteProperty<AngleComponent,GeometryComponent>(constrained_symbol,
                                                            AngleComponent::theta(),
                                                            map),
         p0(point0), p1(point1), p2(point2)
{
    intra_molecule_points = Point::areIntraMoleculePoints(p0,p1) and
                            Point::areIntraMoleculePoints(p0,p2);
}

/** Construct to set the value of 'constrained_symbol' equal to the 
    expression based on the angles within the three points
    'point0', 'point1' and 'point2' */
AngleComponent::AngleComponent(const Symbol &constrained_symbol,
                               const PointRef &point0, const PointRef &point1,
                               const PointRef &point2,
                               const Expression &geometry_expression,
                               const PropertyMap &map)
       : ConcreteProperty<AngleComponent,GeometryComponent>(constrained_symbol,
                                                            geometry_expression,
                                                            map),
         p0(point0), p1(point1), p2(point2)
{
    intra_molecule_points = Point::areIntraMoleculePoints(p0,p1) and
                            Point::areIntraMoleculePoints(p0,p2);
}
  
/** Copy constructor */
AngleComponent::AngleComponent(const AngleComponent &other)
                  : ConcreteProperty<AngleComponent,GeometryComponent>(other),
                    p0(other.p0), p1(other.p1), p2(other.p2),
                    intra_molecule_points(other.intra_molecule_points)
{}

/** Destructor */
AngleComponent::~AngleComponent()
{}

/** Copy assignment operator */
AngleComponent& AngleComponent::operator=(const AngleComponent &other)
{
    if (this != &other)
    {
        GeometryComponent::operator=(other);
        p0 = other.p0;
        p1 = other.p1;
        p2 = other.p2;
        intra_molecule_points = other.intra_molecule_points;
    }
    
    return *this;
}

/** Comparison operator */
bool AngleComponent::operator==(const AngleComponent &other) const
{
    return this == &other or
           (p0 == other.p0 and p1 == other.p1 and p2 == other.p2 and
            GeometryComponent::operator==(other));
}

/** Comparison operator */
bool AngleComponent::operator!=(const AngleComponent &other) const
{
    return not AngleComponent::operator==(other);
}

const char* AngleComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AngleComponent>() );
}

QString AngleComponent::toString() const
{
    return QObject::tr("AngleComponent( %1, %3-%4-%5, %2 )")
                .arg(component().toString(), expression().toString() )
                .arg(p0.read().toString(), p1.read().toString())
                .arg(p2.read().toString());
}

/** Set the space used by the points and distance calculation */
void AngleComponent::setSpace(const Space &new_space)
{
    if (space().equals(new_space) or intra_molecule_points)
        return;

    AngleComponent old_state(*this);
    
    try
    {
        p0.edit().setSpace(new_space);
        p1.edit().setSpace(new_space);
        p2.edit().setSpace(new_space);
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
const Point& AngleComponent::point(int i) const
{
    i = Index(i).map( nPoints() );
    
    switch (i % 3)
    {
        case 0:
            return p0.read();
        case 1:
            return p1.read();
        default:
            return p2.read();
    }
}

/** Return the first point between which the angle is calculated */
const Point& AngleComponent::point0() const
{
    return p0.read();
}

/** Return the second point between which the angle is calculated */
const Point& AngleComponent::point1() const
{
    return p1.read();
}

/** Return the third point between which the angle is calculated */
const Point& AngleComponent::point2() const
{
    return p2.read();
}

/** Return the number of points (3) */
int AngleComponent::nPoints() const
{
    return 3;
}

Triangle AngleComponent::getTriangle() const
{
    if (intra_molecule_points)
    {
        return Triangle(p0.read().point(), p1.read().point(),
                        p2.read().point());
    }
    else
    {
        return Triangle(space().getMinimumImage(p0.read().point(), p1.read().point()),
                        p1.read().point(),
                        space().getMinimumImage(p2.read().point(), p1.read().point()));
    }
}

bool AngleComponent::wouldChange(const Delta &delta, quint32 last_subversion) const
{
    if (delta.hasMoleculeChangeSince(last_subversion))
    {
        return delta.sinceChanged(p0.read(), last_subversion) or
               delta.sinceChanged(p1.read(), last_subversion) or
               delta.sinceChanged(p2.read(), last_subversion);
    }
    else
        return false;
}

Values AngleComponent::getValues(const System &system)
{
    if (p0.read().wouldUpdate(system))
        p0.edit().update(system);

    if (p1.read().wouldUpdate(system))
        p1.edit().update(system);
        
    if (p2.read().wouldUpdate(system))
        p2.edit().update(system);
        
    Triangle trig = getTriangle();

    Values vals;
    vals.set( theta(), trig.angle().to(radians) );
    vals.set( theta012(), trig.angle1().to(radians) );
    vals.set( theta021(), trig.angle2().to(radians) );
    vals.set( theta102(), trig.angle0().to(radians) );
    vals.set( r01(), trig.line2().length() );
    vals.set( r12(), trig.line0().length() );
    vals.set( r02(), trig.line1().length() );
        
    return vals;
}

