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

#ifndef SIREMOL_MOLNUM_H
#define SIREMOL_MOLNUM_H

#include "SireID/number.h"

#include "molid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MolNum;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MolNum&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MolNum&);

namespace SireMol
{

/** This ID number is used to identify a Molecule by
    a unique, program-supplied ID number

    @author Christopher Woods
*/
class SIREMOL_EXPORT MolNum : public SireID::Number, public MolID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const MolNum&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, MolNum&);

public:
    MolNum();
    explicit MolNum(quint32 num);
    
    MolNum(const MolNum &other);

    ~MolNum();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolNum::typeName();
    }
    
    MolNum* clone() const;
    
    static MolNum getUniqueNumber();
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    MolNum& operator=(const MolNum &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const MolNum &other) const;
    bool operator!=(const MolNum &other) const;

    bool operator<(const MolNum &other) const;
    bool operator<=(const MolNum &other) const;
    bool operator>(const MolNum &other) const;
    bool operator>=(const MolNum &other) const;
    
    QList<MolNum> map(const Molecules &molecules) const;
    QList<MolNum> map(const MoleculeGroup&) const;
    QList<MolNum> map(const MolGroupsBase&) const;
};

}

Q_DECLARE_METATYPE(SireMol::MolNum);

SIRE_EXPOSE_CLASS( SireMol::MolNum )

SIRE_END_HEADER

#endif
