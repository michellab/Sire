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

#ifndef SIREMOL_CHAINRESID_H
#define SIREMOL_CHAINRESID_H

#include "chainidentifier.h"
#include "residentifier.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ChainResID;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainResID&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainResID&);

namespace SireMol
{

/** This class holds a combination of a ChainID with a ResID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainResID : public ResID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ChainResID&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ChainResID&);

public:
    ChainResID();
    
    ChainResID(const ChainID &chainid, const ResID &resid);
    
    ChainResID(const ChainResID &other);
    
    ~ChainResID();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ChainResID::typeName();
    }

    ChainResID* clone() const;

    QString toString() const;

    bool operator==(const ChainResID &other) const
    {
        return chainid == other.chainid and resid == other.resid;
    }

    bool operator!=(const ChainResID &other) const
    {
        return not this->operator==(other);
    }

    bool operator==(const SireID::ID &other) const
    {
        return SireID::ID::compare<ChainResID>(*this, other);
    }

    using SireID::ID::operator!=;

    uint hash() const
    {
        return (chainid.hash() << 16) | (resid.hash() & 0x0000FFFF);
    }

    bool isNull() const
    {
        return chainid.isNull() and resid.isNull();
    }

    QList<ResIdx> map(const MolInfo &molinfo) const;

private:
    /** Identifier for the chain */
    ChainIdentifier chainid;
    
    /** Identifier for the residue */
    ResIdentifier resid;
};

}

Q_DECLARE_METATYPE(SireMol::ChainResID);

SIRE_EXPOSE_CLASS( SireMol::ChainResID )

SIRE_END_HEADER

#endif

