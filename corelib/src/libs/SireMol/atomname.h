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

#ifndef SIREMOL_ATOMNAME_H
#define SIREMOL_ATOMNAME_H

#include "SireID/name.h"

#include "atomid.h"
#include "atomidx.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomName;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomName&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomName&);

namespace SireMol
{

/** This class holds the name of an atom. This can be used
    to identify an atom within a residue.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomName : public SireID::Name, public AtomID
{

friend QDataStream& ::operator<<(QDataStream&, const AtomName&);
friend QDataStream& ::operator>>(QDataStream&, AtomName&);

public:
    AtomName();
    
    explicit AtomName(const QString &name);
    
    AtomName(const QString &name, SireID::CaseSensitivity case_sensitivity);
    
    AtomName(const AtomName &other);
    
    ~AtomName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AtomName::typeName();
    }
    
    AtomName* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    AtomName& operator=(const AtomName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const AtomName &other) const;
    
    bool operator!=(const AtomName &other) const;

    QList<AtomIdx> map(const MolInfo &molinfo) const;
};

}

Q_DECLARE_METATYPE(SireMol::AtomName);

SIRE_EXPOSE_CLASS( SireMol::AtomName )

SIRE_END_HEADER

#endif
