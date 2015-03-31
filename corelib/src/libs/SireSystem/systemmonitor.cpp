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

#include <QMutex>

#include "systemmonitor.h"
#include "system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of SystemMonitor
/////////

static const RegisterMetaType<SystemMonitor> r_sysmon( MAGIC_ONLY,
                                                       "SireSystem::SystemMonitor" );
                                                        
/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, 
                                          const SystemMonitor &sysmon)
{
    writeHeader(ds, r_sysmon, 1);
    
    ds << static_cast<const Property&>(sysmon);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, SystemMonitor &sysmon)
{
    VersionID v = readHeader(ds, r_sysmon);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(sysmon);
    }
    else
        throw version_error(v, "1", r_sysmon, CODELOC);
        
    return ds;
}

/** Constructor */
SystemMonitor::SystemMonitor() : Property()
{}

/** Copy constructor */
SystemMonitor::SystemMonitor(const SystemMonitor &other)
           : Property(other)
{}

/** Destructor */
SystemMonitor::~SystemMonitor()
{}

/////////
///////// Implementation of NullMonitor
/////////

static const RegisterMetaType<NullMonitor> r_nullmonitor;
                                                        
/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, 
                                          const NullMonitor &nullmonitor)
{
    writeHeader(ds, r_nullmonitor, 1);
    
    ds << static_cast<const SystemMonitor&>(nullmonitor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, NullMonitor &nullmonitor)
{
    VersionID v = readHeader(ds, r_nullmonitor);
    
    if (v == 1)
    {
        ds >> static_cast<SystemMonitor&>(nullmonitor);
    }
    else
        throw version_error(v, "1", r_nullmonitor, CODELOC);
        
    return ds;
}

/** Constructor */
NullMonitor::NullMonitor() : ConcreteProperty<NullMonitor,SystemMonitor>()
{}

/** Copy constructor */
NullMonitor::NullMonitor(const NullMonitor &other)
            : ConcreteProperty<NullMonitor,SystemMonitor>(other)
{}

/** Destructor */
NullMonitor::~NullMonitor()
{}

/** Copy assignment operator */
NullMonitor& NullMonitor::operator=(const NullMonitor &other)
{
    SystemMonitor::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullMonitor::operator==(const NullMonitor &other) const
{
    return true;
}

/** Comparison operator */
bool NullMonitor::operator!=(const NullMonitor &other) const
{
    return false;
}

/** A null monitor doesn't monitor anything! */
void NullMonitor::monitor(System &system)
{
    return;
}

/** There are no statistics to clear */
void NullMonitor::clearStatistics()
{}

static SharedPolyPointer<NullMonitor> shared_null;

const NullMonitor& SystemMonitor::null()
{
    if (shared_null.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null.constData() == 0)
            shared_null = new NullMonitor();
    }
    
    return *(shared_null.constData());
}

const char* NullMonitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullMonitor>() );
}
