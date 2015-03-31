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

#ifndef SIRESYSTEM_MONITORID_H
#define SIRESYSTEM_MONITORID_H

#include <QList>

#include "SireID/id.h"

#include "SireID/specify.hpp"
#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"

SIRE_BEGIN_HEADER

namespace SireSystem
{

using SireID::Specify;
using SireID::IDAndSet;
using SireID::IDOrSet;

class MonitorIdx;
class MonitorIdentifier;
class MonitorName;

class SystemMonitors;

/** The base class of all system monitor identifiers

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorID : public SireID::ID
{
public:
    typedef MonitorName Index;
    typedef MonitorIdentifier Identifier;
    typedef SystemMonitors SearchObject;

    MonitorID();
    MonitorID(const MonitorID &other);

    virtual ~MonitorID();
    
    static const char* typeName()
    {
        return "SireSystem::MonitorID";
    }
    
    virtual MonitorID* clone() const=0;

    Specify<MonitorID> operator[](int i) const;
    Specify<MonitorID> operator()(int i) const;
    Specify<MonitorID> operator()(int i, int j) const;
    
    IDAndSet<MonitorID> operator+(const MonitorID &other) const;
    IDAndSet<MonitorID> operator&&(const MonitorID &other) const;
    IDAndSet<MonitorID> operator&(const MonitorID &other) const;
    
    IDOrSet<MonitorID> operator*(const MonitorID &other) const;
    IDOrSet<MonitorID> operator||(const MonitorID &other) const;
    IDOrSet<MonitorID> operator|(const MonitorID &other) const;
    
    virtual QList<MonitorName> map(const SystemMonitors &monitors) const=0;

protected:
    QList<MonitorName> processMatches(QList<MonitorName> &matches,
                                      const SystemMonitors &monitors) const;
};

}

#include "monitoridentifier.h"

SIRE_EXPOSE_CLASS( SireSystem::MonitorID )
SIRE_EXPOSE_ALIAS( SireID::Specify<SireSystem::MonitorID>, 
                   SireSystem::Specify_MonitorID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireSystem::MonitorID>, 
                   SireSystem::IDAndSet_MonitorID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireSystem::MonitorID>, 
                   SireSystem::IDOrSet_MonitorID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::Specify<SireSystem::MonitorID>;
template class SireID::IDAndSet<SireSystem::MonitorID>;
template class SireID::IDOrSet<SireSystem::MonitorID>;
#endif

SIRE_END_HEADER

#endif
