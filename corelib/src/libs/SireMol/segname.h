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

#ifndef SIREMOL_SEGNAME_H
#define SIREMOL_SEGNAME_H

#include "SireID/name.h"

#include "segid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class SegName;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SegName&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SegName&);

namespace SireMol
{

/** This class holds the name of a CutGroup.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT SegName : public SireID::Name, public SegID
{

friend QDataStream& ::operator<<(QDataStream&, const SegName&);
friend QDataStream& ::operator>>(QDataStream&, SegName&);

public:
    SegName();
    
    explicit SegName(const QString &name);
    
    SegName(const QString &name, SireID::CaseSensitivity case_sensitivity);
    
    SegName(const SegName &other);
    
    ~SegName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return SegName::typeName();
    }
    
    SegName* clone() const
    {
        return new SegName(*this);
    }
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    SegName& operator=(const SegName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const SegName &other) const;
    
    bool operator!=(const SegName &other) const;

    QList<SegIdx> map(const MolInfo &molinfo) const;
};

}

Q_DECLARE_METATYPE(SireMol::SegName);

SIRE_EXPOSE_CLASS( SireMol::SegName )

SIRE_END_HEADER

#endif
