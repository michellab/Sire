/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_RESNUM_H
#define SIREMOL_RESNUM_H

#include "SireID/number.h"

#include "resid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ResNum;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResNum&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResNum&);

namespace SireMol
{

/** This ID number is used to identify a CutGroup by the user-supplied
    number

    @author Christopher Woods
*/
class SIREMOL_EXPORT ResNum : public SireID::Number, public ResID
{

friend QDataStream& ::operator<<(QDataStream&, const ResNum&);
friend QDataStream& ::operator>>(QDataStream&, ResNum&);

public:
    ResNum();

    explicit ResNum(quint32 num);

    ResNum(const ResNum &other);

    ~ResNum();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ResNum::typeName();
    }
    
    ResNum* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    ResNum& operator=(const ResNum &other);
    
    bool operator==(const SireID::ID &other) const;
    bool operator==(const ResNum &other) const;
    bool operator!=(const ResNum &other) const;

    bool operator<(const ResNum &other) const;
    bool operator<=(const ResNum &other) const;
    bool operator>(const ResNum &other) const;
    bool operator>=(const ResNum &other) const;

    QList<ResIdx> map(const MolInfo &molinfo) const;
};

}

Q_DECLARE_METATYPE(SireMol::ResNum);

SIRE_EXPOSE_CLASS( SireMol::ResNum )

SIRE_END_HEADER

#endif
