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

#ifndef SIRESYSTEM_SYSNAME_H
#define SIRESYSTEM_SYSNAME_H

#include "SireID/name.h"

#include "sysid.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class SysName;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::SysName&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::SysName&);

namespace SireSystem
{

/** This class holds the name of a simulation system
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT SysName : public SireID::Name, public SysID
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const SysName&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, SysName&);

public:
    SysName();
    explicit SysName(const QString &name);
    
    SysName(const SysName &other);
    
    ~SysName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return SysName::typeName();
    }
    
    SysName* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    SysName& operator=(const SysName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const SysName &other) const;
    
    bool operator!=(const SysName &other) const;
    
    QList<SysIdx> map(const Systems &systems) const;
};

}

Q_DECLARE_METATYPE( SireSystem::SysName );

SIRE_EXPOSE_CLASS( SireSystem::SysName )

SIRE_END_HEADER

#endif
