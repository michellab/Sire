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

#ifndef SIREMOL_CHAINNAME_H
#define SIREMOL_CHAINNAME_H

#include "SireID/name.h"

#include "chainid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ChainName;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainName&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainName&);

namespace SireMol
{

/** This class holds the name of an atom. This can be used
    to identify an atom within a residue.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainName : public SireID::Name, public ChainID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ChainName&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ChainName&);

public:
    ChainName();
    
    explicit ChainName(const QString &name);
    
    ChainName(const QString &name, SireID::CaseSensitivity case_sensitivity);
    
    ChainName(const ChainName &other);
    
    ~ChainName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ChainName::typeName();
    }
    
    ChainName* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    ChainName& operator=(const ChainName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const ChainName &other) const;
    
    bool operator!=(const ChainName &other) const;

    QList<ChainIdx> map(const MolInfo &molinfo) const;
};

}

Q_DECLARE_METATYPE(SireMol::ChainName);

SIRE_EXPOSE_CLASS( SireMol::ChainName )

SIRE_END_HEADER

#endif

