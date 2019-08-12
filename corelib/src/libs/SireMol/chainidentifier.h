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

#ifndef SIREMOL_CHAINIDENTIFIER_H
#define SIREMOL_CHAINIDENTIFIER_H

#include "chainid.h"

#include <boost/shared_ptr.hpp>

namespace SireMol
{
class ChainIdentifier;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainIdentifier&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainIdentifier&);

namespace SireMol
{

class ChainIdx;

/** This is the polymorphic holder of all Chain IDs */
class SIREMOL_EXPORT ChainIdentifier : public ChainID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ChainIdentifier&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ChainIdentifier&);

public:
    ChainIdentifier();
    ChainIdentifier(const ChainID &chainid);
    ChainIdentifier(const ChainIdentifier &other);
    
    ~ChainIdentifier();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ChainIdentifier::typeName();
    }
    
    ChainIdentifier* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const ChainID& base() const;
    
    ChainIdentifier& operator=(const ChainIdentifier &other);
    ChainIdentifier& operator=(const ChainID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const ChainIdentifier &other) const;
    bool operator!=(const ChainIdentifier &other) const;
    
    bool operator==(const ChainID &other) const;
    bool operator!=(const ChainID &other) const;
    
    QList<ChainIdx> map(const MolInfo &molinfo) const;

private:
    /** Pointer to the ChainID */
    boost::shared_ptr<ChainID> d;
};

SIRE_ALWAYS_INLINE uint qHash(const ChainIdentifier &chainid)
{
    return chainid.hash();
}

}

#include "chainidx.h"
#include "residx.h"
#include "atomidx.h"

Q_DECLARE_METATYPE( SireID::Specify<SireMol::ChainID> );
Q_DECLARE_METATYPE( SireMol::AtomsIn<SireMol::ChainID> );
Q_DECLARE_METATYPE( SireMol::ResIn<SireMol::ChainID> );
Q_DECLARE_METATYPE( SireID::IDAndSet<SireMol::ChainID> );
Q_DECLARE_METATYPE( SireID::IDOrSet<SireMol::ChainID> );
Q_DECLARE_METATYPE( SireID::MatchAll<SireMol::ChainID> );
Q_DECLARE_METATYPE( SireID::InvertMatch<SireMol::ChainID> );

Q_DECLARE_METATYPE(SireMol::ChainIdentifier);

#endif
