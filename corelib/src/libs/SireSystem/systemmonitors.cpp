/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "systemmonitors.h"

#include "system.h"

#include "monitorname.h"
#include "monitoridx.h"

#include "tostring.h"

#include "SireSystem/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireStream;

static const RegisterMetaType<SystemMonitors> r_sysmons(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SystemMonitors &sysmons)
{
    writeHeader(ds, r_sysmons, 1);
    
    SharedDataStream sds(ds);

    sds << sysmons.mons_by_name << sysmons.mons_by_idx
        << sysmons.mons_by_frequency << sysmons.stepnum;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SystemMonitors &sysmons)
{
    VersionID v = readHeader(ds, r_sysmons);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> sysmons.mons_by_name >> sysmons.mons_by_idx
            >> sysmons.mons_by_frequency >> sysmons.stepnum;
    }
    else
        throw version_error(v, "1", r_sysmons, CODELOC);

    return ds;
}

/** Constructor */
SystemMonitors::SystemMonitors() : stepnum(0)
{}

/** Copy constructor */
SystemMonitors::SystemMonitors(const SystemMonitors &other)
               : mons_by_name(other.mons_by_name),
                 mons_by_idx(other.mons_by_idx),
                 mons_by_frequency(other.mons_by_frequency),
                 stepnum(other.stepnum)
{}

/** Destructor */
SystemMonitors::~SystemMonitors()
{}

/** Copy assignment operator */
SystemMonitors& SystemMonitors::operator=(const SystemMonitors &other)
{
    if (this != &other)
    {
        mons_by_name = other.mons_by_name;
        mons_by_idx = other.mons_by_idx;
        mons_by_frequency = other.mons_by_frequency;
        stepnum = other.stepnum;
    }
    
    return *this;
}

/** Comparison operator */
bool SystemMonitors::operator==(const SystemMonitors &other) const
{
    return mons_by_name == other.mons_by_name and
           mons_by_idx == other.mons_by_idx and
           mons_by_frequency == other.mons_by_frequency and
           stepnum == other.stepnum;
}

/** Comparison operator */
bool SystemMonitors::operator!=(const SystemMonitors &other) const
{
    return not this->operator==(other);
}

/** Return whether or not this is empty (contains no monitors) */
bool SystemMonitors::isEmpty() const
{
    return mons_by_idx.isEmpty();
}
    
/** Return the name of the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
MonitorName SystemMonitors::monitorName(const MonitorID &monid) const
{
    QList<MonitorName> names = monid.map(*this);
    
    if (names.count() > 1)
        throw SireSystem::duplicate_monitor( QObject::tr(
            "More than one system monitor matches the ID '%1'. Matching "
            "monitors have names %2.")
                .arg(monid.toString(), Sire::toString(names)),
                    CODELOC );

    return names.at(0);
}
      
/** Return the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
const SystemMonitor& SystemMonitors::operator[](const MonitorID &monid) const
{
    return mons_by_name.constFind( this->monitorName(monid) )->read();
}

/** Return the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
const SystemMonitor& SystemMonitors::at(const MonitorID &monid) const
{
    return this->operator[](monid);
}

/** Simple shortcut function

    \throw SireSystem::missing_monitor
*/
QList<MonitorName> SystemMonitors::map(const MonitorName &monname) const
{
    if (not mons_by_name.contains(monname))
        throw SireSystem::missing_monitor( QObject::tr(
            "There is no system monitor called %1. Available "
            "monitors are %2.")
                .arg(monname, Sire::toString(mons_by_name.keys())),
                    CODELOC );
                    
    QList<MonitorName> names;
    names.append(monname);
    return names;
}

/** Return the name of the monitor at index 'monidx'

    \throw SireError::invalid_index
*/
QList<MonitorName> SystemMonitors::map(const MonitorIdx &monidx) const
{
    QList<MonitorName> names;
    
    names.append( mons_by_idx.at( monidx.map(mons_by_idx.count()) ) );
    
    return names;
}

/** Return the names of the monitors that match the ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
QList<MonitorName> SystemMonitors::map(const MonitorID &monid) const
{
    return monid.map(*this);
}

/** Return the frequency of the monitor with ID 'monid'

    \throw SireSystem::duplicate_monitor
    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
int SystemMonitors::getFrequency(const MonitorID &monid) const
{
    MonitorName name = this->monitorName(monid);
    
    for (QHash< quint32, QList<MonitorName> >::const_iterator
                                        it = mons_by_frequency.constBegin();
         it != mons_by_frequency.constEnd();
         ++it)
    {
        if (it.value().contains(name))
            return it.key();
    }
    
    return 0;
}

/** Return the list of all monitor names */
QList<MonitorName> SystemMonitors::monitorNames() const
{
    return mons_by_name.keys();
}

/** Return the list of all monitor names */
QList<MonitorName> SystemMonitors::names() const
{
    return this->monitorNames();
}

/** Add a system monitor 'monitor', identified by the name 'name', which
    will be updated every 'frequency' steps.
    
    \throw SireSystem::duplicate_monitor
*/
void SystemMonitors::add(const QString &name, const SystemMonitor &monitor,
                         int frequency)
{
    MonitorName monname( name );

    if (mons_by_name.contains(monname))
        throw SireSystem::duplicate_monitor( QObject::tr(
            "Cannot add the monitor of type %1 as a monitor with "
            "the same name (%2, of type %3) is already present in this "
            "collection.")
                .arg(monitor.what()).arg(name).arg(mons_by_name.value(monname)->what()),
                    CODELOC );
                    
    mons_by_name.insert(monname, monitor);
    mons_by_idx.append(monname);

    if (frequency < 0)
        frequency = 0;
        
    mons_by_frequency[frequency].append(monname);
}

/** Add the monitors from 'other' to this set

    \throw SireSystem::duplicate_monitor
*/
void SystemMonitors::add(const SystemMonitors &other)
{
    if (this->isEmpty())
    {
        this->operator=(other);
        return;
    }
    else if (other.isEmpty())
    {
        return;
    }

    SystemMonitors old_state(*this);
    
    try
    {
        foreach (const MonitorName &name, other.names())
        {
            this->add( name, other[name], other.getFrequency(name) );
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add all of the monitors in 'other' to this set, adding them
    with the frequency 'frequency' 
    
    \throw SireSystem::duplicate_monitor
*/
void SystemMonitors::add(const SystemMonitors &other, int frequency)
{
    SystemMonitors new_monitors(other);
    new_monitors.setAllFrequency(frequency);
    
    this->add(new_monitors);
}

/** Set the frequency of all of the monitors to 'frequency' */
void SystemMonitors::setAllFrequency(int frequency)
{
    //do we need to make any change?
    if (this->isEmpty())
        return;
    
    if (frequency < 0)
        frequency = 0;
                
    if (mons_by_frequency.count() == 1 and 
        mons_by_frequency.contains( quint32(frequency) ))
    {
        return;
    }
    
    mons_by_frequency.clear();
    mons_by_frequency.insert( quint32(frequency), mons_by_idx );
}

/** Remove all of the monitors that match the ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
void SystemMonitors::remove(const MonitorID &monid)
{
    QList<MonitorName> names = monid.map(*this);
    
    foreach (const MonitorName &name, names)
    {
        mons_by_name.remove(name);
        mons_by_idx.removeAll(name);
        
        QMutableHashIterator< quint32,QList<MonitorName> > it(mons_by_frequency);
        
        while (it.hasNext())
        {
            it.next();
            it.value().removeAll(name);
            
            if (it.value().isEmpty())
                it.remove();
        }
    }
}

/** Remove all of the monitors from this set */
void SystemMonitors::removeAll()
{
    mons_by_name.clear();
    mons_by_idx.clear();
    mons_by_frequency.clear();
}

/** Completely clear any statistics held in these monitors */
void SystemMonitors::clearStatistics()
{
    for (QHash<MonitorName,SysMonPtr>::iterator it = mons_by_name.begin();
         it != mons_by_name.end();
         ++it)
    {
        it.value().edit().clearStatistics();
    }
}

/** Completely clear the statistics held by the monitors that
    match the ID 'monid' - this does nothing if there are no
    monitors that match this ID */
void SystemMonitors::clearStatistics(const MonitorID &monid)
{
    QList<MonitorName> monitor_names;
    
    try
    {
        monitor_names = monid.map(*this);
    }
    catch(...)
    {}
    
    foreach (MonitorName monitor_name, monitor_names)
    {
        BOOST_ASSERT( mons_by_name.contains(monitor_name) );
        
        mons_by_name[monitor_name].edit().clearStatistics();
    }
}

/** Set the frequency of all monitors that match the ID 'monid' so that
    they are updated every 'frequency' steps
    
    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
void SystemMonitors::setFrequency(const MonitorID &monid, int frequency)
{
    QList<MonitorName> names = monid.map(*this);
    
    foreach (const MonitorName &name, names)
    {
        QMutableHashIterator< quint32,QList<MonitorName> > it(mons_by_frequency);
        
        while (it.hasNext())
        {
            it.next();
            it.value().removeAll(name);
            
            if (it.value().isEmpty())
                it.remove();
        }
    }
    
    if (frequency < 0)
        frequency = 0;
    
    mons_by_frequency[frequency] += names;
}

/** Return the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
const SystemMonitor& SystemMonitors::monitor(const MonitorID &monid) const
{
    return this->operator[](monid);
}

/** Return all of the monitors that match the ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
QList<SysMonPtr> SystemMonitors::monitors(const MonitorID &monid) const
{
    QList<MonitorName> names = monid.map(*this);
    
    QList<SysMonPtr> mons;
    
    foreach (const MonitorName &name, names)
    {
        mons.append( mons_by_name.value(name) );
    }
    
    return mons;
}

/** Return the list of all monitors in this set, in the order they
    appear in this set */
QList<SysMonPtr> SystemMonitors::monitors() const
{
    QList<SysMonPtr> mons;
    
    foreach (const MonitorName &name, mons_by_idx)
    {
        mons.append( mons_by_name.value(name) );
    }
    
    return mons;
}

/** Return the list of all monitors in this set, in the order they
    appear in this set */
QList<SysMonPtr> SystemMonitors::list() const
{
    return this->monitors();
}

/** Return the number of monitors in this set */
int SystemMonitors::nMonitors() const
{
    return mons_by_name.count();
}

/** Return the number of monitors in this set */
int SystemMonitors::count() const
{
    return this->nMonitors();
}

/** Return the number of monitors in this set */
int SystemMonitors::size() const
{
    return this->nMonitors();
}

/** Update the monitors by monitoring the system 'system' */
void SystemMonitors::monitor(System &system)
{
    SystemMonitors old_state(*this);

    try
    {

        //increment the step number
        stepnum += 1;
    
        for (QHash< quint32,QList<MonitorName> >::const_iterator 
                                            it = mons_by_frequency.constBegin();
             it != mons_by_frequency.constEnd();
             ++it)
        {
            if (stepnum % it.key() == 0)
            {
                //it is time to update these monitors!
                foreach (const MonitorName &name, it.value())
                {
                    mons_by_name[name].edit().monitor(system);
                }
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

const char* SystemMonitors::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SystemMonitors>() );
}
