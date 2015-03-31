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

#include "geometrycomponent.h"
#include "delta.h"

#include "SireVol/cartesian.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireCAS;
using namespace SireVol;
using namespace SireStream;

static const RegisterMetaType<GeometryComponent> r_geomcomp( MAGIC_ONLY,
                                                      GeometryComponent::typeName() );
                                                      
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds,
                                          const GeometryComponent &geomcomp)
{
    writeHeader(ds, r_geomcomp, 1);
    
    SharedDataStream sds(ds);

    PropertyMap map;
    map["space"] = geomcomp.space_property;
    
    sds << geomcomp.constrained_symbol << geomcomp.geometry_expression
        << map
        << static_cast<const Constraint&>(geomcomp);
        
    return ds;
}

QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, GeometryComponent &geomcomp)
{
    VersionID v = readHeader(ds, r_geomcomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        geomcomp.clearLastSystem();
        
        PropertyMap map;
        
        sds >> geomcomp.constrained_symbol >> geomcomp.geometry_expression
            >> map >> static_cast<Constraint&>(geomcomp);
            
        geomcomp.expected_value = 0;
        geomcomp.current_values = Values();
        geomcomp.spce = Cartesian();
        geomcomp.space_property = map["space"];
    }
    else
        throw version_error(v, "1", r_geomcomp, CODELOC);
        
    return ds;
}

/** Null constructor */
GeometryComponent::GeometryComponent(const PropertyMap &map) 
                  : Constraint(), expected_value(0), space_property( map["space"] )
{}

/** Internal constructor */
GeometryComponent::GeometryComponent(const Symbol &symbol, const Expression &expression,
                                     const PropertyMap &map)
                  : Constraint(), constrained_symbol(symbol),
                    expected_value(0), geometry_expression(expression),
                    space_property( map["space"] )
{}

/** Copy constructor */
GeometryComponent::GeometryComponent(const GeometryComponent &other)
                  : Constraint(other),
                    constrained_symbol(other.constrained_symbol),
                    expected_value(other.expected_value),
                    current_values(other.current_values),
                    geometry_expression(other.geometry_expression),
                    space_property(other.space_property), spce(other.spce)
{}

/** Destructor */
GeometryComponent::~GeometryComponent()
{}

/** Copy assignment operator */
GeometryComponent& GeometryComponent::operator=(const GeometryComponent &other)
{
    if (this != &other)
    {
        Constraint::operator=(other);
        constrained_symbol = other.constrained_symbol;
        expected_value = other.expected_value;
        current_values = other.current_values;
        geometry_expression = other.geometry_expression;
        space_property = other.space_property;
        spce = other.spce;
    }
    
    return *this;
}

/** Comparison operator */
bool GeometryComponent::operator==(const GeometryComponent &other) const
{
    return this == &other or
           (Constraint::operator==(other) and
            constrained_symbol == other.constrained_symbol and 
            geometry_expression == other.geometry_expression and
            space_property == other.space_property);
}

/** Comparison operator */
bool GeometryComponent::operator!=(const GeometryComponent &other) const
{
    return not GeometryComponent::operator==(other);
}

const char* GeometryComponent::typeName()
{
    return "SireSystem::GeometryComponent";
}

/** Return the symbol representing the constrained component */
const Symbol& GeometryComponent::component() const
{
    return constrained_symbol;
}

/** Return the expression used to calculate the constrained value
    from the geometry */
const Expression& GeometryComponent::expression() const
{
    return geometry_expression;
}

/** Internal function used to set the space used to evaluate the geometry */
void GeometryComponent::setSpace(const Space &space)
{
    spce = space;
}

const Space& GeometryComponent::space() const
{
    return spce.read();
}

void GeometryComponent::setSystem(const System &system)
{
    if (Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system))
        return;

    if (system.containsProperty(space_property))
        this->setSpace( system.property(space_property).asA<Space>() );
    else
        this->setSpace( Cartesian() );
    
    Values vals = this->getValues(system);
    expected_value = geometry_expression(vals);
    
    if (system.hasConstantComponent(constrained_symbol))
    {
        Constraint::setSatisfied( system, system.constant(constrained_symbol) == 
                                                                     expected_value );
    }
    else
        Constraint::setSatisfied( system, false );
}

bool GeometryComponent::mayChange(const Delta &delta, quint32 last_subversion) const
{
    return delta.sinceChanged(space_property, last_subversion) or
           delta.sinceChanged(constrained_symbol, last_subversion) or
           this->wouldChange(delta, last_subversion);
}

bool GeometryComponent::fullApply(Delta &delta)
{
    this->setSystem(delta.deltaSystem());
    
    if (not Constraint::wasLastSatisfied())
    {
        return delta.update(constrained_symbol, expected_value);
    }
    else
    {
        return false;
    }
}

bool GeometryComponent::deltaApply(Delta &delta, quint32 last_subversion)
{
    if ( delta.sinceChanged(space_property, last_subversion) )
    {
        return this->fullApply(delta);
    }
    else if (this->wouldChange(delta, last_subversion))
    {
        Values vals = this->getValues( delta.deltaSystem() );
        expected_value = geometry_expression(vals);
        
        return delta.update( constrained_symbol, expected_value );
    }
    else if ( delta.sinceChanged(constrained_symbol, last_subversion) )
    {
        return delta.update( constrained_symbol, expected_value );
    }
    else
        return false;
}
