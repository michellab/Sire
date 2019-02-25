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

#ifndef SIRESYSTEM_SYSTEMMONITOR_H
#define SIRESYSTEM_SYSTEMMONITOR_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class SystemMonitor;
class NullMonitor;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::SystemMonitor&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::SystemMonitor&);

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::NullMonitor&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::NullMonitor&);

namespace SireSystem
{

class System;

/** This is the virtual base class of all system monitors. A system
    monitor is an object that monitors a system during a simulation,
    e.g. collecting the average energy, saving a radial distribution
    function etc.
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT SystemMonitor : public SireBase::Property
{
public:
    SystemMonitor();
    
    SystemMonitor(const SystemMonitor &other);
    
    virtual ~SystemMonitor();
    
    static const char* typeName()
    {
        return "SireSystem::SystemMonitor";
    }
    
    virtual SystemMonitor* clone() const=0;
    
    virtual void clearStatistics()=0;
    
    virtual void monitor(System &system)=0;
    
    static const NullMonitor& null();
};

/** This is a null monitor that doesn't monitor anything */
class SIRESYSTEM_EXPORT NullMonitor
           : public SireBase::ConcreteProperty<NullMonitor,SystemMonitor>
{
public:
    NullMonitor();
    
    NullMonitor(const NullMonitor &other);
    
    ~NullMonitor();
    
    NullMonitor& operator=(const NullMonitor &other);
    
    bool operator==(const NullMonitor &other) const;
    bool operator!=(const NullMonitor &other) const;
    
    static const char *typeName();
    
    void clearStatistics();
    
    void monitor(System &system);
};

typedef SireBase::PropPtr<SystemMonitor> SysMonPtr;

}

Q_DECLARE_METATYPE( SireSystem::NullMonitor )

SIRE_EXPOSE_CLASS( SireSystem::SystemMonitor )
SIRE_EXPOSE_CLASS( SireSystem::NullMonitor )

SIRE_EXPOSE_PROPERTY( SireSystem::SysMonPtr, SireSystem::SystemMonitor )

SIRE_END_HEADER

#endif
