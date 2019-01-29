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

#ifndef SIREMOL_MOLNAME_H
#define SIREMOL_MOLNAME_H

#include "SireID/name.h"

#include "molid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MolName;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MolName&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MolName&);

namespace SireMol
{

/** This class holds the name of a Molecule.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MolName : public SireID::Name, public MolID
{

friend QDataStream& ::operator<<(QDataStream&, const MolName&);
friend QDataStream& ::operator>>(QDataStream&, MolName&);

public:
    MolName();
    explicit MolName(const QString &name);
    
    MolName(const QString &name, SireID::CaseSensitivity case_sensitivity);
    
    MolName(const MolName &other);
    
    ~MolName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolName::typeName();
    }
    
    MolName* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    MolName& operator=(const MolName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const MolName &other) const;
    
    bool operator!=(const MolName &other) const;
    
    QList<MolNum> map(const Molecules &molecules) const;
    QList<MolNum> map(const MoleculeGroup &molgroup) const;
    QList<MolNum> map(const MolGroupsBase &molgroups) const;
};

}

Q_DECLARE_METATYPE(SireMol::MolName);

SIRE_EXPOSE_CLASS( SireMol::MolName )

SIRE_END_HEADER

#endif

