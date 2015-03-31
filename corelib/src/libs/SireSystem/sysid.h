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

#ifndef SIRESYSTEM_SYSID_H
#define SIRESYSTEM_SYSID_H

#include <QList>

#include "SireID/id.h"

#include "SireID/specify.hpp"
#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"

SIRE_BEGIN_HEADER

namespace SireSystem
{

class SysID;
class SysIdx;
class SysIdentifier;
class SysName;

using SireID::Specify;
using SireID::IDAndSet;
using SireID::IDOrSet;

/** Dummy set */
class SIRESYSTEM_EXPORT Systems
{
public:
    Systems()
    {}
    
    ~Systems()
    {}
    
    QList<SysIdx> getSystems() const;
    
    QList<SysIdx> map(const SysID &sysid) const;
};

/** The base class of all system identifiers

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT SysID : public SireID::ID
{
public:
    typedef SysIdx Index;
    typedef SysIdentifier Identifier;
    typedef Systems SearchObject;

    SysID();
    SysID(const SysID &other);

    virtual ~SysID();
    
    static const char* typeName()
    {
        return "SireSystem::SysID";
    }
    
    virtual SysID* clone() const=0;
    
    Specify<SysID> operator[](int i) const;
    Specify<SysID> operator()(int i) const;
    Specify<SysID> operator()(int i, int j) const;
    
    IDAndSet<SysID> operator+(const SysID &other) const;
    IDAndSet<SysID> operator&&(const SysID &other) const;
    IDAndSet<SysID> operator&(const SysID &other) const;
    
    IDOrSet<SysID> operator*(const SysID &other) const;
    IDOrSet<SysID> operator||(const SysID &other) const;
    IDOrSet<SysID> operator|(const SysID &other) const;
    
    virtual QList<SysIdx> map(const Systems &systems) const=0;

protected:
    QList<SysIdx> processMatches(QList<SysIdx> &matches, 
                                 const Systems &systems) const;
};

}

#include "sysidentifier.h"

SIRE_EXPOSE_CLASS( SireSystem::SysID )
SIRE_EXPOSE_ALIAS( SireID::Specify<SireSystem::SysID>, SireSystem::Specify_SysID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireSystem::SysID>, SireSystem::IDAndSet_SysID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireSystem::SysID>, SireSystem::IDOrSet_SysID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::Specify<SireSystem::SysID>;
template class SireID::IDAndSet<SireSystem::SysID>;
template class SireID::IDOrSet<SireSystem::SysID>;
#endif

SIRE_END_HEADER

#endif
