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

#include "monitormonitor.h"

#include "SireSystem/system.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<MonitorMonitor> r_monmon;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                          const MonitorMonitor &monmon)
{
    writeHeader(ds, r_monmon, 1);
    
    SharedDataStream sds(ds);
    
    sds << monmon.monitor_states << monmon.monitor_id
        << monmon.clear_original << monmon.remove_original
        << static_cast<const SystemMonitor&>(monmon);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          MonitorMonitor &monmon)
{
    VersionID v = readHeader(ds, r_monmon);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> monmon.monitor_states >> monmon.monitor_id
            >> monmon.clear_original >> monmon.remove_original
            >> static_cast<SystemMonitor&>(monmon);
    }
    else
        throw version_error(v, "1", r_monmon, CODELOC);
        
    return ds;
}

/** Null constructor */
MonitorMonitor::MonitorMonitor() 
               : ConcreteProperty<MonitorMonitor,SystemMonitor>(),
                 clear_original(false), remove_original(false)
{}
               
/** Construct to monitor the SystemMonitor identified by 'id', 
    optionally clearing the original monitor if 'clear_original' is
    true, and optionally removing the original monitor
    if 'remove_original' is true */
MonitorMonitor::MonitorMonitor(const MonitorID &id, bool clear, bool remove)
               : ConcreteProperty<MonitorMonitor,SystemMonitor>(),
                 monitor_id(id), clear_original(clear),
                 remove_original(remove)
{}

/** Copy constructor */
MonitorMonitor::MonitorMonitor(const MonitorMonitor &other)
               : ConcreteProperty<MonitorMonitor,SystemMonitor>(other),
                 monitor_states(other.monitor_states),
                 monitor_id(other.monitor_id), clear_original(other.clear_original),
                 remove_original(other.remove_original)
{}

/** Destructor */
MonitorMonitor::~MonitorMonitor()
{}

/** Copy assignment operator */
MonitorMonitor& MonitorMonitor::operator=(const MonitorMonitor &other)
{
    if (this != &other)
    {
        SystemMonitor::operator=(other);
        
        monitor_states = other.monitor_states;
        monitor_id = other.monitor_id;
        clear_original = other.clear_original;
        remove_original = other.remove_original;
    }
    
    return *this;
}

/** Comparison operator */
bool MonitorMonitor::operator==(const MonitorMonitor &other) const
{
    return (this == &other) or
           (monitor_states == other.monitor_states and
            monitor_id == other.monitor_id and
            clear_original == other.clear_original and
            remove_original == other.remove_original and
            SystemMonitor::operator==(other));
}

/** Comparison operator */
bool MonitorMonitor::operator!=(const MonitorMonitor &other) const
{
    return not this->operator==(other);
}

/** Return the number of states monitored so far */
int MonitorMonitor::nStates() const
{
    return monitor_states.count();
}

/** Return the number of states monitored so far */
int MonitorMonitor::count() const
{
    return this->nStates();
}

/** Return the number of states monitored so far */
int MonitorMonitor::size() const
{
    return this->nStates();
}

/** Return the ith state monitored

    \throw SireError::invalid_index
*/
const SystemMonitor& MonitorMonitor::operator[](int i) const
{
    return monitor_states.at( Index(i).map(monitor_states.count()) );
}

/** Return the ith state monitored

    \throw SireError::invalid_index
*/
const SystemMonitor& MonitorMonitor::at(int i) const
{
    return this->operator[](i);
}

/** Return all of the states monitored */
const QList<SysMonPtr>& MonitorMonitor::states() const
{
    return monitor_states;
}

/** Set whether or not to clear the statistics of the original
    monitor when this monitor takes the copy */
void MonitorMonitor::setClearOriginal(bool clear)
{
    clear_original = clear;
}

/** Set whether or not to remove the original monitor when
    this monitor takes a copy (effectively thus moving 
    the monitor from the system to this MonitorMonitor) */
void MonitorMonitor::setRemoveOriginal(bool remove)
{
    remove_original = remove;
}

/** Return whether or not this MonitorMonitor will clear the
    original monitor whenever it takes a copy */
bool MonitorMonitor::clearOriginal() const
{
    return clear_original;
}

/** Return whether or not this MonitorMonitor will remove
    the original monitor from the system whenever it
    takes a copy */
bool MonitorMonitor::removeOriginal() const
{
    return remove_original;
}

/** Completely clear statistics */
void MonitorMonitor::clearStatistics()
{
    monitor_states.clear();
}

/** Monitor the passed system */
void MonitorMonitor::monitor(System &system)
{
    monitor_states.append( system.monitor(monitor_id) );
    
    if (remove_original)
        system.remove(monitor_id);
        
    else if (clear_original)
        system.clearStatistics(monitor_id);
}

const char* MonitorMonitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MonitorMonitor>() );
}
