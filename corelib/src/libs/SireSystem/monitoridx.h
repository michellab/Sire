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

#ifndef SIRESYSTEM_MONITORIDX_H
#define SIRESYSTEM_MONITORIDX_H

#include "SireID/index.h"

#include "monitorid.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MonitorIdx;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::MonitorIdx&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::MonitorIdx&);

namespace SireSystem
{

/** This is an ID object that is used to index system monitors (e.g. index
    in a list or array).

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorIdx 
                : public SireID::Index_T_<MonitorIdx>, public MonitorID
{

friend QDataStream& ::operator<<(QDataStream&, const MonitorIdx&);
friend QDataStream& ::operator>>(QDataStream&, MonitorIdx&);

public:
    MonitorIdx();
    explicit MonitorIdx(qint32 idx);
    
    MonitorIdx(const MonitorIdx &other);
    
    ~MonitorIdx();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MonitorIdx::typeName();
    }
    
    MonitorIdx* clone() const;
    
    static MonitorIdx null();
    
    bool isNull() const;
    
    uint hash() const;

    QString toString() const;
    
    MonitorIdx& operator=(const MonitorIdx &other);
    
    bool operator==(const SireID::ID &other) const;
    
    using SireID::Index_T_<MonitorIdx>::operator=;

    using SireID::Index_T_<MonitorIdx>::operator==;
    using SireID::Index_T_<MonitorIdx>::operator!=;

    using SireID::Index_T_<MonitorIdx>::operator+=;
    using SireID::Index_T_<MonitorIdx>::operator++;
    using SireID::Index_T_<MonitorIdx>::operator-=;
    using SireID::Index_T_<MonitorIdx>::operator--;
    
    using SireID::Index_T_<MonitorIdx>::map;
    
    QList<MonitorName> map(const SystemMonitors &monitors) const;
};
    
}

Q_DECLARE_TYPEINFO(SireSystem::MonitorIdx, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(SireSystem::MonitorIdx);

SIRE_EXPOSE_CLASS( SireSystem::MonitorIdx )

SIRE_END_HEADER

#endif
