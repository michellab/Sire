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

#ifndef SIRESYSTEM_MONITORS_H
#define SIRESYSTEM_MONITORS_H

#include "systemmonitor.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class SystemMonitors;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::SystemMonitors&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::SystemMonitors&);

namespace SireSystem
{

class MonitorID;
class MonitorIdx;
class MonitorName;

/** This class holds a set of SystemMonitor objects, and controls
    when those monitors are evaluated on a system
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT SystemMonitors
{

friend QDataStream& ::operator<<(QDataStream&, const SystemMonitors&);
friend QDataStream& ::operator>>(QDataStream&, SystemMonitors&);

public:
    SystemMonitors();
    
    SystemMonitors(const SystemMonitors &other);
    
    ~SystemMonitors();

    static const char* typeName();
    
    const char* what() const
    {
        return SystemMonitors::typeName();
    }

    SystemMonitors& operator=(const SystemMonitors &other);
    
    bool operator==(const SystemMonitors &other) const;
    bool operator!=(const SystemMonitors &other) const;
    
    const SystemMonitor& operator[](const MonitorID &monid) const;

    const SystemMonitor& at(const MonitorID &monid) const;
    
    bool isEmpty() const;

    QList<MonitorName> map(const MonitorName &monname) const;
    QList<MonitorName> map(const MonitorIdx &monidx) const;
    QList<MonitorName> map(const MonitorID &monid) const;

    MonitorName monitorName(const MonitorID &monid) const;

    QList<MonitorName> monitorNames() const;

    QList<MonitorName> names() const;
    QList<SysMonPtr> list() const;

    void add(const QString &name, const SystemMonitor &monitor,
             int frequency = 1);
    
    void add(const SystemMonitors &other);
    void add(const SystemMonitors &other, int frequency);
    
    void remove(const MonitorID &monid);
    void removeAll();
    
    void clearStatistics();
    void clearStatistics(const MonitorID &monid);
    
    void setAllFrequency(int frequency);
    void setFrequency(const MonitorID &monid, int frequency);
   
    int getFrequency(const MonitorID &monid) const;
     
    const SystemMonitor& monitor(const MonitorID &monid) const;
    
    QList<SysMonPtr> monitors(const MonitorID &monid) const;
    
    QList<SysMonPtr> monitors() const;
    
    int nMonitors() const;
    int count() const;
    int size() const;
    
    void monitor(System &system);
    
private:
    /** All of the monitors, indexed by name */
    QHash<MonitorName,SysMonPtr> mons_by_name;
    
    /** The names of all of the monitors in the order they
        appear in this collection */
    QList<MonitorName> mons_by_idx;
    
    /** The frequency of each monitor */
    QHash< quint32, QList<MonitorName> > mons_by_frequency;
    
    /** The current step number of this set of monitors */
    quint32 stepnum;
};

}

Q_DECLARE_METATYPE( SireSystem::SystemMonitors )

SIRE_EXPOSE_CLASS( SireSystem::SystemMonitors )

SIRE_END_HEADER

#endif
