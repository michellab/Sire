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

#ifndef SIREMOL_MGNUM_H
#define SIREMOL_MGNUM_H

#include "SireID/number.h"

#include "mgid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MGNum;
}

QDataStream& operator<<(QDataStream&, const SireMol::MGNum&);
QDataStream& operator>>(QDataStream&, SireMol::MGNum&);

namespace SireMol
{

/** This ID number is used to identify a MoleculeGroup by
    a unique, program-supplied ID number

    @author Christopher Woods
*/
class SIREMOL_EXPORT MGNum : public SireID::Number, public MGID
{

friend QDataStream& ::operator<<(QDataStream&, const MGNum&);
friend QDataStream& ::operator>>(QDataStream&, MGNum&);

public:
    MGNum();
    explicit MGNum(quint32 num);
    
    MGNum(const MGNum &other);

    ~MGNum();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MGNum::typeName();
    }
    
    MGNum* clone() const;
    
    static MGNum getUniqueNumber();
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    MGNum& operator=(const MGNum &other);
    
    bool operator==(const SireID::ID &other) const;
    bool operator==(const MGNum &other) const;
    bool operator!=(const MGNum &other) const;

    bool operator<(const MGNum &other) const;
    bool operator<=(const MGNum &other) const;
    bool operator>(const MGNum &other) const;
    bool operator>=(const MGNum &other) const;
    
    QList<MGNum> map(const MolGroupsBase&) const;
};

}

Q_DECLARE_METATYPE(SireMol::MGNum);

SIRE_EXPOSE_CLASS( SireMol::MGNum )

SIRE_END_HEADER

#endif
