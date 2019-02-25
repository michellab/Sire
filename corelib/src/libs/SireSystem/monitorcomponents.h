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

#ifndef SIRESYSTEM_MONITORCOMPONENTS_H
#define SIRESYSTEM_MONITORCOMPONENTS_H

#include "systemmonitor.h"

#include "SireMaths/accumulator.h"

#include "SireCAS/symbol.h"
#include "SireCAS/symbols.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MonitorComponents;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::MonitorComponents&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::MonitorComponents&);

namespace SireSystem
{

class MonitorComponent;

using SireMaths::Accumulator;
using SireMaths::Average;

using SireCAS::Symbol;
using SireCAS::Symbols;

/** This is a monitor that can be used to monitor large numbers 
    of components of the system (of even all components). This
    monitor is similar to MonitorComponent, but is better suited
    to situations where you want to monitor everything (or 
    everything except for a small number of components), where
    it would be messy to specify lots of individual MonitorComponent
    monitors
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorComponents
            : public SireBase::ConcreteProperty<MonitorComponents,SystemMonitor>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const MonitorComponents&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, MonitorComponents&);

public:
    MonitorComponents();
    
    MonitorComponents(const Accumulator &accumulator);
    
    MonitorComponents(const Symbol &component, 
                      const Accumulator &accumulator = Average());
    MonitorComponents(const Symbols &components, 
                      const Accumulator &accumulator = Average());
    
    MonitorComponents(const MonitorComponent &component_monitor);
    
    MonitorComponents(const MonitorComponents &other);
    
    ~MonitorComponents();
    
    MonitorComponents& operator=(const MonitorComponents &other);
    
    bool operator==(const MonitorComponents &other) const;
    bool operator!=(const MonitorComponents &other) const;
    
    static const char* typeName();

    void excludeComponent(const Symbol &component);
    void excludeComponent(const Symbols &components);

    const Symbols& includeComponents() const;
    const Symbols& excludeComponents() const;
    
    Symbols monitoredComponents() const;
    
    const Accumulator& accumulatorTemplate() const;
    
    const Accumulator& accumulator(const Symbol &component) const;
    
    void clearStatistics();
    
    void monitor(System &system);
    
private:
    /** The symbols representing components to be monitored.
        If this is empty, then all components will be monitored */
    Symbols include_symbols;
    
    /** The symbols representing components that must
        not be monitored */
    Symbols exclude_symbols;
    
    /** The accumulator used to accumulate components - this
        is a template */
    SireMaths::AccumulatorPtr accumulator_template;
    
    /** The actual accumulators associated with each component 
        that has been monitored */
    QHash<Symbol,SireMaths::AccumulatorPtr> accumes;
};

}

Q_DECLARE_METATYPE( SireSystem::MonitorComponents )

SIRE_EXPOSE_CLASS( SireSystem::MonitorComponents )

SIRE_END_HEADER

#endif
