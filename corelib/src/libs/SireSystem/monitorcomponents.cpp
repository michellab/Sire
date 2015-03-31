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

#include "monitorcomponents.h"

#include "monitorcomponent.h"
#include "system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireCAS/errors.h"

#include "tostring.h"

using namespace SireSystem;
using namespace SireCAS;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<MonitorComponents> r_moncomps;

/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds,
                                          const MonitorComponents &moncomps)
{
    writeHeader(ds, r_moncomps, 1);
    
    SharedDataStream sds(ds);
    
    sds << moncomps.include_symbols << moncomps.exclude_symbols
        << moncomps.accumulator_template << moncomps.accumes
        << static_cast<const SystemMonitor&>(moncomps);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, 
                                          MonitorComponents &moncomps)
{
    VersionID v = readHeader(ds, r_moncomps);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> moncomps.include_symbols >> moncomps.exclude_symbols
            >> moncomps.accumulator_template >> moncomps.accumes
            >> static_cast<SystemMonitor&>(moncomps);
    }
    else
        throw version_error(v, "1", r_moncomps, CODELOC);
        
    return ds;
}

/** Null constructor */
MonitorComponents::MonitorComponents()
                  : ConcreteProperty<MonitorComponents,SystemMonitor>()
{}

/** Construct a monitor to monitor all components using the
    accumulator 'accumulator' */
MonitorComponents::MonitorComponents(const Accumulator &accumulator)
                  : ConcreteProperty<MonitorComponents,SystemMonitor>(),
                    accumulator_template(accumulator)
{}

/** Construct a monitor to monitor the component 'component' using
    the accumulator 'accumulator' */
MonitorComponents::MonitorComponents(const Symbol &component, 
                                     const Accumulator &accumulator)
                  : ConcreteProperty<MonitorComponents,SystemMonitor>(),
                    include_symbols(component),
                    accumulator_template(accumulator)
{}

/** Construct a monitor to monitor the components 'components'
    using the accumulator 'accumulator' */
MonitorComponents::MonitorComponents(const Symbols &components, 
                                     const Accumulator &accumulator)
                  : ConcreteProperty<MonitorComponents,SystemMonitor>(),
                    include_symbols(components),
                    accumulator_template(accumulator)
{}

/** Construct a MonitorComponents that is the (near) equivalent of
    the MonitorComponent 'component_monitor' */
MonitorComponents::MonitorComponents(const MonitorComponent &component_monitor)
                  : ConcreteProperty<MonitorComponents,SystemMonitor>(),
                    include_symbols(component_monitor.component()),
                    accumulator_template(component_monitor.accumulator())
{}

/** Copy constructor */
MonitorComponents::MonitorComponents(const MonitorComponents &other)
                  : ConcreteProperty<MonitorComponents,SystemMonitor>(other),
                    include_symbols(other.include_symbols),
                    exclude_symbols(other.exclude_symbols),
                    accumulator_template(other.accumulator_template),
                    accumes(other.accumes)
{}

/** Destructor */
MonitorComponents::~MonitorComponents()
{}

/** Copy assignment operator */
MonitorComponents& MonitorComponents::operator=(const MonitorComponents &other)
{
    if (this != &other)
    {
        include_symbols = other.include_symbols;
        exclude_symbols = other.exclude_symbols;
        accumulator_template = other.accumulator_template;
        accumes = other.accumes;
    }
    
    return *this;
}

/** Comparison operator */
bool MonitorComponents::operator==(const MonitorComponents &other) const
{
    return (this == &other) or
           (include_symbols == other.include_symbols and 
            exclude_symbols == other.exclude_symbols and 
            accumulator_template == other.accumulator_template and
            accumes == other.accumes and
            SystemMonitor::operator==(other) );
}

/** Comparison operator */
bool MonitorComponents::operator!=(const MonitorComponents &other) const
{
    return not this->operator==(other);
}

/** Make sure that the component 'component' is not monitored. */
void MonitorComponents::excludeComponent(const Symbol &component)
{
    exclude_symbols.insert(component);
}

/** Make sure that the components in 'components' are not monitored */
void MonitorComponents::excludeComponent(const Symbols &components)
{
    exclude_symbols.insert(components);
}

/** Return the components that will be monitored
    (if they exist, and not if they are excluded). If this
    is empty, then all components are monitored */
const Symbols& MonitorComponents::includeComponents() const
{
    return include_symbols;
}

/** Return the components that will definitely not be monitored */
const Symbols& MonitorComponents::excludeComponents() const
{
    return exclude_symbols;
}

/** Return the set of symbols that have been monitored so far
    (so have valid accumulators) */
Symbols MonitorComponents::monitoredComponents() const
{
    return accumes.keys();
}

/** Return the accumulator that is the template used for new accumulators
    that are created when a new component is monitored */
const Accumulator& MonitorComponents::accumulatorTemplate() const
{
    return accumulator_template;
}

/** Return the accumulator for the component 'component'

    \throw SireCAS::missing_symbol
*/
const Accumulator& MonitorComponents::accumulator(const Symbol &component) const
{
    QHash<Symbol,AccumulatorPtr>::const_iterator it = accumes.constFind(component);
    
    if (it == accumes.constEnd())
        throw SireCAS::missing_symbol( QObject::tr(
            "There is no accumulator for the component represented by the symbol %1. "
            "Available components are %2.")
                .arg(component.toString())
                .arg(Sire::toString(accumes.keys())), CODELOC );

    return it.value();
}

/** Completely clear the statistics */
void MonitorComponents::clearStatistics()
{
    accumes.clear();
}

/** Monitor the system 'system' - this will only accumulate symbols
    that represent existing components - this does nothing for components
    that don't exist in the system */
void MonitorComponents::monitor(System &system)
{
    if (accumulator_template.isNull())
        //there is no accumulator
        return;

    Values vals;

    if (exclude_symbols.isEmpty())
    {
        if (include_symbols.isEmpty())
        {
            vals = system.energies();
        }
        else
        {
            Symbols available_components = include_symbols;
            available_components.intersect( system.componentSymbols() );
            
            if (not available_components.isEmpty())
                vals = system.componentValues(available_components);
        }
    }
    else
    {
        Symbols available_components = system.componentSymbols();
        
        if (not include_symbols.isEmpty())
            available_components.intersect(include_symbols);
            
        available_components.subtract(exclude_symbols);
        
        if (not available_components.isEmpty())
            vals = system.energies(available_components);
    }
    
    const QHash<SymbolID,double> &values = vals.values();
    
    for (QHash<SymbolID,double>::const_iterator it = values.constBegin();
         it != values.constEnd();
         ++it)
    {
        Symbol component(it.key());
        
        if (not accumes.contains(component))
            accumes.insert(component, accumulator_template);
            
        accumes[component].edit().accumulate(it.value());
    }
}

const char* MonitorComponents::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MonitorComponents>() );
}
