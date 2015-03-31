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

#ifndef SIREFF_FFIDX_H
#define SIREFF_FFIDX_H

#include "SireID/index.h"

#include "ffid.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class FFIdx;
}

QDataStream& operator<<(QDataStream&, const SireFF::FFIdx&);
QDataStream& operator>>(QDataStream&, SireFF::FFIdx&);

namespace SireFF
{

/** This is an ID object that is used to index forcefields (e.g. index
    in a list or array).

    @author Christopher Woods
*/
class SIREFF_EXPORT FFIdx : public SireID::Index_T_<FFIdx>, public FFID
{

friend QDataStream& ::operator<<(QDataStream&, const FFIdx&);
friend QDataStream& ::operator>>(QDataStream&, FFIdx&);

public:
    FFIdx();
    explicit FFIdx(qint32 idx);
    
    FFIdx(const FFIdx &other);
    
    ~FFIdx();
    
    static const char* typeName();
    
    const char* what() const
    {
        return FFIdx::typeName();
    }
    
    FFIdx* clone() const;
    
    static FFIdx null();
    
    bool isNull() const;
    
    uint hash() const;

    QString toString() const;
    
    FFIdx& operator=(const FFIdx &other);
    
    bool operator==(const SireID::ID &other) const;
    
    using SireID::Index_T_<FFIdx>::operator=;

    using SireID::Index_T_<FFIdx>::operator==;
    using SireID::Index_T_<FFIdx>::operator!=;

    using SireID::Index_T_<FFIdx>::operator+=;
    using SireID::Index_T_<FFIdx>::operator++;
    using SireID::Index_T_<FFIdx>::operator-=;
    using SireID::Index_T_<FFIdx>::operator--;
    
    using SireID::Index_T_<FFIdx>::map;
    
    QList<FFIdx> map(const ForceFields &ffields) const;
};
    
}

Q_DECLARE_METATYPE(SireFF::FFIdx);

SIRE_EXPOSE_CLASS( SireFF::FFIdx )

SIRE_END_HEADER

#endif
