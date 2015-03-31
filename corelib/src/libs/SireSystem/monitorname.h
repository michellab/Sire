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

#ifndef SIRESYSTEM_MONITORNAME_H
#define SIRESYSTEM_MONITORNAME_H

#include "SireID/name.h"

#include "monitorid.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MonitorName;
}

QDataStream& operator<<(QDataStream&, const SireSystem::MonitorName&);
QDataStream& operator>>(QDataStream&, SireSystem::MonitorName&);

namespace SireSystem
{

/** This class holds the name of a simulation system
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorName : public SireID::Name, public MonitorID
{

friend QDataStream& ::operator<<(QDataStream&, const MonitorName&);
friend QDataStream& ::operator>>(QDataStream&, MonitorName&);

public:
    MonitorName();
    explicit MonitorName(const QString &name);
    
    MonitorName(const MonitorName &other);
    
    ~MonitorName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MonitorName::typeName();
    }
    
    MonitorName* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    MonitorName& operator=(const MonitorName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const MonitorName &other) const;
    
    bool operator!=(const MonitorName &other) const;
    
    QList<MonitorName> map(const SystemMonitors &monitors) const;
};

}

Q_DECLARE_METATYPE( SireSystem::MonitorName );

SIRE_EXPOSE_CLASS( SireSystem::MonitorName )

SIRE_END_HEADER

#endif
