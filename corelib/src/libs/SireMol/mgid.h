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

#ifndef SIREMOL_MGID_H
#define SIREMOL_MGID_H

#include "SireID/id.h"

#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"
#include "SireID/specify.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{

using SireID::IDAndSet;
using SireID::IDOrSet;
using SireID::Specify;

class MGIdx;
class MGIdentifier;
class MGNum;

class MolGroupsBase;

/** This is the base class of all identifiers that are used 
    to identify a MoleculeGroup

    @author Christopher Woods
*/
class SIREMOL_EXPORT MGID : public SireID::ID
{
public:
    typedef MGNum Index;
    typedef MGIdentifier Identifier;
    typedef MolGroupsBase SearchObject;

    MGID();
    
    MGID(const MGID &other);
    
    virtual ~MGID();

    static const char* typeName()
    {
        return "SireMol::MGID";
    }
    
    Specify<MGID> operator[](int i) const;
    Specify<MGID> operator()(int i) const;
    Specify<MGID> operator()(int i, int j) const;
    
    IDAndSet<MGID> operator+(const MGID &other) const;
    IDOrSet<MGID> operator*(const MGID &other) const;

    IDAndSet<MGID> operator&&(const MGID &other) const;
    IDAndSet<MGID> operator&(const MGID &other) const;

    IDOrSet<MGID> operator||(const MGID &other) const;
    IDOrSet<MGID> operator|(const MGID &other) const;

    virtual MGID* clone() const=0;

    virtual QList<MGNum> map(const MolGroupsBase &molgroups) const=0;

protected:
    void processMatches(QList<MGNum> &matches, const MolGroupsBase &molgroups) const;
};

}

#include "mgidentifier.h"

SIRE_EXPOSE_CLASS( SireMol::MGID )
SIRE_EXPOSE_ALIAS( SireID::Specify<SireMol::MGID>, SireMol::Specify_MGID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireMol::MGID>, SireMol::IDAndSet_MGID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireMol::MGID>, SireMol::IDOrSet_MGID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::IDAndSet<SireMol::MGID>;
template class SireID::IDOrSet<SireMol::MGID>;
template class SireID::Specify<SireMol::MGID>;
#endif

SIRE_END_HEADER

#endif
