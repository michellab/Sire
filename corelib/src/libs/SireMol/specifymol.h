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

#ifndef SIREMOL_SPECIFYMOL_H
#define SIREMOL_SPECIFYMOL_H

#include "molid.h"
#include "molidentifier.h"

#include "SireID/index.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class SpecifyMol;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SpecifyMol&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SpecifyMol&);

namespace SireMol
{

/** This class allow for the specification of specific
    matching molecules
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT SpecifyMol : public MolID
{

friend QDataStream& ::operator<<(QDataStream&, const SpecifyMol&);
friend QDataStream& ::operator>>(QDataStream&, SpecifyMol&);

public:
    SpecifyMol();
    SpecifyMol(const MolID &molid);
    SpecifyMol(const MolID &molid, int i);
    SpecifyMol(const MolID &molid, int i, int j);
    
    SpecifyMol(const SpecifyMol &other);
    
    ~SpecifyMol();
    
    static const char* typeName();
    
    const char* what() const
    {
        return SpecifyMol::typeName();
    }
    
    SpecifyMol* clone() const;
    
    SpecifyMol& operator=(const SpecifyMol &other);
    
    bool operator==(const SpecifyMol &other) const;
    bool operator!=(const SpecifyMol &other) const;
    
    bool operator==(const SireID::ID &other) const;
    bool operator!=(const SireID::ID &other) const;
    
    uint hash() const;
    
    QString toString() const;
    
    bool isNull() const;
    
    QList<MolNum> map(const Molecules &molecules) const;
    QList<MolNum> map(const MoleculeGroup &molgroup) const;
    QList<MolNum> map(const MolGroupsBase &molgroups) const;

private:
    /** The underlying molecule ID */
    MolIdentifier molid;
    
    /** The range of molecules to match */
    SireID::Index strt, end;
};

}

Q_DECLARE_METATYPE( SireMol::SpecifyMol )

SIRE_EXPOSE_CLASS( SireMol::SpecifyMol )

SIRE_END_HEADER

#endif
