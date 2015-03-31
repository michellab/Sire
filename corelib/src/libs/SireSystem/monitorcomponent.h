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

#ifndef SIRESYSTEM_MONITORCOMPONENT_H
#define SIRESYSTEM_MONITORCOMPONENT_H

#include "systemmonitor.h"

#include "SireMaths/accumulator.h"
#include "SireCAS/symbol.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MonitorComponent;
}

QDataStream& operator<<(QDataStream&, const SireSystem::MonitorComponent&);
QDataStream& operator>>(QDataStream&, SireSystem::MonitorComponent&);

namespace SireSystem 
{

using SireMaths::Accumulator;

using SireCAS::Symbol;

/** This monitor is used to monitor the value of a component of the system.
    It can be used to generate an average of that component. It could be
    used, for example, to generate average energies, or accumulate
    free energies during a simulation
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorComponent 
        : public SireBase::ConcreteProperty<MonitorComponent,SystemMonitor>
{

friend QDataStream& ::operator<<(QDataStream&, const MonitorComponent&);
friend QDataStream& ::operator>>(QDataStream&, MonitorComponent&);

public:
    MonitorComponent();
    MonitorComponent(const Symbol &component);
    MonitorComponent(const Symbol &component, 
                     const Accumulator &accumulator);
                     
    MonitorComponent(const MonitorComponent &other);
    
    ~MonitorComponent();
    
    MonitorComponent& operator=(const MonitorComponent &other);
    
    static const char* typeName();
    
    bool operator==(const MonitorComponent &other) const;
    bool operator!=(const MonitorComponent &other) const;
    
    const Symbol& component() const;
    const Accumulator& accumulator() const;

    void clearStatistics();

    void monitor(System &system);

private:
    /** The symbol representing the component being monitored */
    Symbol monitored_component;
    
    /** The accumulator that is used to accumulate the average
        (and other properties) */
    SireMaths::AccumulatorPtr accume;
};

}

Q_DECLARE_METATYPE( SireSystem::MonitorComponent )

SIRE_EXPOSE_CLASS( SireSystem::MonitorComponent )

SIRE_END_HEADER

#endif
