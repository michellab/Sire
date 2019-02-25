/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "monitorcomponent.h"
#include "system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireCAS;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<MonitorComponent> r_moncomp;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const MonitorComponent &moncomp)
{
    writeHeader(ds, r_moncomp, 1);
    
    SharedDataStream sds(ds);
    
    sds << moncomp.monitored_component
        << moncomp.accume
        << static_cast<const SystemMonitor&>(moncomp);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MonitorComponent &moncomp)
{
    VersionID v = readHeader(ds, r_moncomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> moncomp.monitored_component
            >> moncomp.accume
            >> static_cast<SystemMonitor&>(moncomp);
    }
    else
        throw version_error( v, "1", r_moncomp, CODELOC );
        
    return ds;
}

/** Null constructor */
MonitorComponent::MonitorComponent()
                 : ConcreteProperty<MonitorComponent,SystemMonitor>()
{}

/** Construct a monitor to collect the average value of the component
    represented by the symbol 'component' */
MonitorComponent::MonitorComponent(const Symbol &component)
                 : ConcreteProperty<MonitorComponent,SystemMonitor>(),
                   monitored_component(component),
                   accume( Average() )
{}

/** Construct a monitor to accumulate the value of the component
    represented by the symbol 'component' using the accumulator
    in 'accumulator' */
MonitorComponent::MonitorComponent(const Symbol &component, 
                                   const Accumulator &accumulator)
                 : ConcreteProperty<MonitorComponent,SystemMonitor>(),
                   monitored_component(component),
                   accume(accumulator)
{}

/** Copy constructor */                 
MonitorComponent::MonitorComponent(const MonitorComponent &other)
                 : ConcreteProperty<MonitorComponent,SystemMonitor>(other),
                   monitored_component(other.monitored_component),
                   accume(other.accume)
{}

/** Destructor */
MonitorComponent::~MonitorComponent()
{}

/** Copy assignment operator */
MonitorComponent& MonitorComponent::operator=(const MonitorComponent &other)
{
    if (this != &other)
    {
        monitored_component = other.monitored_component;
        accume = other.accume;
        SystemMonitor::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool MonitorComponent::operator==(const MonitorComponent &other) const
{
    return monitored_component == other.monitored_component and
           accume == other.accume;
}

/** Comparison operator */
bool MonitorComponent::operator!=(const MonitorComponent &other) const
{
    return not this->operator==(other);
}

/** Return the symbol representing the component being monitored */
const Symbol& MonitorComponent::component() const
{
    return monitored_component;
}

/** Return the accumulator that is being used to accumulate the
    values of the component being monitored */
const Accumulator& MonitorComponent::accumulator() const
{
    return accume;
}

/** Clear the statistics in this monitor */
void MonitorComponent::clearStatistics()
{
    accume.edit().clear();
}

/** Call this function to add the statistics of the monitored
    component from the passed system to the accumulator */
void MonitorComponent::monitor(System &system)
{
    accume.edit().accumulate( system.componentValue(monitored_component) );
}

const char* MonitorComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MonitorComponent>() );
}
