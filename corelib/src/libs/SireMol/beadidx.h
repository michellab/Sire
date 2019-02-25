/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMOL_BEADIDX_H
#define SIREMOL_BEADIDX_H

#include "SireID/index.h"

#include "beadid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class BeadIdx;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::BeadIdx&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::BeadIdx&);

namespace SireMol
{

class CGAtomIdx;

/** This is an ID object that is used to index CutGroups

    @author Christopher Woods
*/
class SIREMOL_EXPORT BeadIdx 
       : public SireID::Index_T_<BeadIdx>, public BeadID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const BeadIdx&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, BeadIdx&);

public:
    BeadIdx();
    
    explicit BeadIdx(qint32 idx);
    
    BeadIdx(const BeadIdx &other);
    
    ~BeadIdx();
    
    static const char* typeName();

    const char* what() const
    {
        return SireID::Index_T_<BeadIdx>::what();
    }
    
    BeadIdx* clone() const;
    
    static BeadIdx null();
    
    bool isNull() const;
    
    uint hash() const;

    QString toString() const;
    
    BeadIdx& operator=(const BeadIdx &other);
    
    bool operator==(const SireID::ID &other) const;
    
    using SireID::Index_T_<BeadIdx>::operator=;

    using SireID::Index_T_<BeadIdx>::operator==;
    using SireID::Index_T_<BeadIdx>::operator!=;

    using SireID::Index_T_<BeadIdx>::operator+=;
    using SireID::Index_T_<BeadIdx>::operator++;
    using SireID::Index_T_<BeadIdx>::operator-=;
    using SireID::Index_T_<BeadIdx>::operator--;
    
    using SireID::Index_T_<BeadIdx>::map;
};
    
}

Q_DECLARE_METATYPE(SireMol::BeadIdx);

SIRE_EXPOSE_CLASS( SireMol::BeadIdx )

SIRE_END_HEADER

#endif
