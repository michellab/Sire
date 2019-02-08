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

#ifndef SIREMOL_RESIDX_H
#define SIREMOL_RESIDX_H

#include "SireID/index.h"

#include "resid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ResIdx;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResIdx&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResIdx&);

namespace SireMol
{

/** This is an ID object that is used to index CutGroups

    @author Christopher Woods
*/
class SIREMOL_EXPORT ResIdx 
       : public SireID::Index_T_<ResIdx>, public ResID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ResIdx&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ResIdx&);

public:
    ResIdx();
    
    explicit ResIdx(quint32 idx);
    
    ResIdx(const ResIdx &other);
    
    ~ResIdx();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ResIdx::typeName();
    }
    
    ResIdx* clone() const;
    
    static ResIdx null();
    
    bool isNull() const;
    
    uint hash() const;

    QString toString() const;
    
    ResIdx& operator=(const ResIdx &other);
    
    bool operator==(const SireID::ID &other) const;
    
    using SireID::Index_T_<ResIdx>::operator=;

    using SireID::Index_T_<ResIdx>::operator==;
    using SireID::Index_T_<ResIdx>::operator!=;

    using SireID::Index_T_<ResIdx>::operator+=;
    using SireID::Index_T_<ResIdx>::operator++;
    using SireID::Index_T_<ResIdx>::operator-=;
    using SireID::Index_T_<ResIdx>::operator--;
    
    using SireID::Index_T_<ResIdx>::map;
    
    QList<ResIdx> map(const MolInfo &molinfo) const;
};
    
}

Q_DECLARE_METATYPE(SireMol::ResIdx);

SIRE_EXPOSE_CLASS( SireMol::ResIdx )

SIRE_END_HEADER

#endif

