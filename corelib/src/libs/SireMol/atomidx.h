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

#ifndef SIREMOL_ATOMIDX_H
#define SIREMOL_ATOMIDX_H

#include "SireID/index.h"

#include "atomid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomIdx;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomIdx&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomIdx&);

namespace SireMol
{

class MoleculeInfo;

/** This is an ID object that is used to index atoms (e.g. index
    in a list or array, or in a molecule).

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomIdx : public SireID::Index_T_<AtomIdx>, public AtomID
{

friend QDataStream& ::operator<<(QDataStream&, const AtomIdx&);
friend QDataStream& ::operator>>(QDataStream&, AtomIdx&);

public:
    AtomIdx();
    explicit AtomIdx(qint32 idx);
    
    AtomIdx(const AtomIdx &other);
    
    ~AtomIdx();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AtomIdx::typeName();
    }
    
    AtomIdx* clone() const;
    
    static AtomIdx null();
    
    bool isNull() const;
    
    uint hash() const;

    QString toString() const;
    
    AtomIdx& operator=(const AtomIdx &other);
    
    bool operator==(const SireID::ID &other) const;
    
    using SireID::Index_T_<AtomIdx>::operator=;

    using SireID::Index_T_<AtomIdx>::operator==;
    using SireID::Index_T_<AtomIdx>::operator!=;

    using SireID::Index_T_<AtomIdx>::operator+=;
    using SireID::Index_T_<AtomIdx>::operator++;
    using SireID::Index_T_<AtomIdx>::operator-=;
    using SireID::Index_T_<AtomIdx>::operator--;
    
    using SireID::Index_T_<AtomIdx>::map;
    
    QList<AtomIdx> map(const MolInfo &molinfo) const;
};
    
}

Q_DECLARE_METATYPE(SireMol::AtomIdx);

SIRE_EXPOSE_CLASS( SireMol::AtomIdx )

SIRE_END_HEADER

#endif
