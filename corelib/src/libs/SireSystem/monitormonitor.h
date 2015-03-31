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

#ifndef SIRESYSTEM_MONITOR_MONITOR_H
#define SIRESYSTEM_MONITOR_MONITOR_H

#include "systemmonitor.h"
#include "monitoridentifier.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MonitorMonitor;
}

QDataStream& operator<<(QDataStream&, const SireSystem::MonitorMonitor&);
QDataStream& operator>>(QDataStream&, SireSystem::MonitorMonitor&);

namespace SireSystem
{

/** This is a monitor that can be used to monitor other monitors.
    It is useful to use to save the current state of a monitor,
    e.g. snap-shotting the average during the simulation, or
    snap-shotting the PDB. Equally, it can be used to transfer
    monitor data between systems, or from a system to a supra-system.
    
    The MonitorMonitor can be non-destructive (the original monitor
    is just copied), or destructive (the original monitor is copied,
    then cleared), or annihilative (the original monitor is copied,
    then completely removed from the system)
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorMonitor
          : public SireBase::ConcreteProperty<MonitorMonitor,SystemMonitor>
{

friend QDataStream& ::operator<<(QDataStream&, const MonitorMonitor&);
friend QDataStream& ::operator>>(QDataStream&, MonitorMonitor&);

public:
    MonitorMonitor();
    MonitorMonitor(const MonitorID &id, bool clear_original=false,
                                        bool remove_original=false);
    
    MonitorMonitor(const MonitorMonitor &other);
    
    ~MonitorMonitor();
    
    MonitorMonitor& operator=(const MonitorMonitor &other);
    
    bool operator==(const MonitorMonitor &other) const;
    bool operator!=(const MonitorMonitor &other) const;
    
    static const char* typeName();

    const SystemMonitor& operator[](int i) const;
    
    const SystemMonitor& at(int i) const;

    int size() const;
    int count() const;
    int nStates() const;

    const QList<SysMonPtr>& states() const;

    void setClearOriginal(bool clear);
    void setRemoveOriginal(bool remove);
    
    bool clearOriginal() const;
    bool removeOriginal() const;

    void clearStatistics();
    
    void monitor(System &system);
   
private:
    /** The monitor being monitored */
    QList<SysMonPtr> monitor_states;
    
    /** The ID of the monitor to monitor */
    MonitorIdentifier monitor_id;

    /** Whether or not to clear the original monitor */
    bool clear_original;
    
    /** Whether or not to remove the original monitor */
    bool remove_original;
};

}

Q_DECLARE_METATYPE( SireSystem::MonitorMonitor )

SIRE_EXPOSE_CLASS( SireSystem::MonitorMonitor )

SIRE_END_HEADER

#endif
