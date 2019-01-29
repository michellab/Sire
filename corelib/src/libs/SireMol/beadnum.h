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

#ifndef SIREMOL_BEADNUM_H
#define SIREMOL_BEADNUM_H

#include "SireID/number.h"

#include "beadid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class BeadNum;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::BeadNum&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::BeadNum&);

namespace SireMol
{

/** This ID number is used to identify a Bead by the user-supplied
    number

    @author Christopher Woods
*/
class SIREMOL_EXPORT BeadNum : public SireID::Number, public BeadID
{

friend QDataStream& ::operator<<(QDataStream&, const BeadNum&);
friend QDataStream& ::operator>>(QDataStream&, BeadNum&);

public:
    BeadNum();

    explicit BeadNum(quint32 num);

    BeadNum(const BeadNum &other);

    ~BeadNum();
    
    static const char* typeName();
    
    const char* what() const
    {
        return BeadNum::typeName();
    }
    
    BeadNum* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    BeadNum& operator=(const BeadNum &other);
    
    bool operator==(const SireID::ID &other) const;
    bool operator==(const BeadNum &other) const;
    bool operator!=(const BeadNum &other) const;

    bool operator<(const BeadNum &other) const;
    bool operator<=(const BeadNum &other) const;
    bool operator>(const BeadNum &other) const;
    bool operator>=(const BeadNum &other) const;
};

}

Q_DECLARE_METATYPE(SireMol::BeadNum);

SIRE_EXPOSE_CLASS( SireMol::BeadNum )

SIRE_END_HEADER

#endif
