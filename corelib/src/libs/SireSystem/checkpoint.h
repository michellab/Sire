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

#ifndef SIRESYSTEM_CHECKPOINT_H
#define SIRESYSTEM_CHECKPOINT_H

#include "system.h"

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class CheckPoint;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::CheckPoint&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::CheckPoint&);

namespace SireSystem
{

/** This class holds a checkpoint of a system. This allows you to 
    save the current state of a system so that you can restore
    it at a later point (or even save it to disk/database)
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT CheckPoint
          : public SireBase::ConcreteProperty<CheckPoint,SireBase::Property>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const CheckPoint&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, CheckPoint&);

public:
    CheckPoint();
    CheckPoint(const System &system);
    
    CheckPoint(const CheckPoint &other);
    
    ~CheckPoint();
    
    static const char* typeName();
    
    CheckPoint& operator=(const System &system);
    CheckPoint& operator=(const CheckPoint &other);
    
    bool operator==(const CheckPoint &other) const;
    bool operator!=(const CheckPoint &other) const;
    
    operator System() const;

private:
    /** A copy of the system that has been checkpointed */
    System old_system;
};

}

Q_DECLARE_METATYPE( SireSystem::CheckPoint )

SIRE_EXPOSE_CLASS( SireSystem::CheckPoint )

SIRE_END_HEADER

#endif
