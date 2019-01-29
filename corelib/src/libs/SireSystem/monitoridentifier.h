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

#ifndef SIRESYSTEM_MONITORIDENTIFIER_H
#define SIRESYSTEM_MONITORIDENTIFIER_H

#include "monitorid.h"

#include <boost/shared_ptr.hpp>

namespace SireSystem
{
class MonitorIdentifier;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::MonitorIdentifier&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::MonitorIdentifier&);

namespace SireSystem
{

/** This is a generic holder for any MonitorID class! 

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorIdentifier : public MonitorID
{

friend QDataStream& ::operator<<(QDataStream&, const MonitorIdentifier&);
friend QDataStream& ::operator>>(QDataStream&, MonitorIdentifier&);

public:
    MonitorIdentifier();
    MonitorIdentifier(const MonitorID &monid);
    MonitorIdentifier(const MonitorIdentifier &other);
    
    ~MonitorIdentifier();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MonitorIdentifier::typeName();
    }
    
    bool isNull() const;
    
    MonitorIdentifier* clone() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const MonitorID& base() const;
    
    MonitorIdentifier& operator=(const MonitorIdentifier &other);
    MonitorIdentifier& operator=(const MonitorID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const MonitorIdentifier &other) const;
    bool operator!=(const MonitorIdentifier &other) const;
    
    bool operator==(const MonitorID &other) const;
    bool operator!=(const MonitorID &other) const;

    QList<MonitorName> map(const SystemMonitors &monitors) const;

private:
    /** Pointer to the MonitorID */
    boost::shared_ptr<MonitorID> d;
};

inline uint qHash(const MonitorIdentifier &monid)
{
    return monid.hash();
}

}

#include "monitorname.h"

Q_DECLARE_METATYPE( SireID::Specify<SireSystem::MonitorID> )
Q_DECLARE_METATYPE( SireID::IDAndSet<SireSystem::MonitorID> )
Q_DECLARE_METATYPE( SireID::IDOrSet<SireSystem::MonitorID> )

Q_DECLARE_METATYPE(SireSystem::MonitorIdentifier);

#endif
