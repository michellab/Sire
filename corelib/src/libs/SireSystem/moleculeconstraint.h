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

#ifndef SIRESYSTEM_MOLECULECONSTRAINT_H
#define SIRESYSTEM_MOLECULECONSTRAINT_H

#include <QUuid>

#include "constraint.h"

#include "SireBase/majorminorversion.h"

#include "SireMol/molecules.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MoleculeConstraint;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::MoleculeConstraint&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::MoleculeConstraint&);

namespace SireMol
{
class MolNum;
class Molecules;
}

namespace SireSystem
{

using SireBase::Version;

using SireMol::MolNum;
using SireMol::Molecules;

/** This is the base class of all Molecule constraints. These
    are constraints that affect molecules in a system. A molecule
    constraint works by being updated by a system and returning
    the molecules that must then be changed to maintain the 
    constraint
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MoleculeConstraint : public Constraint
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeConstraint&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeConstraint&);

public:
    MoleculeConstraint();
    MoleculeConstraint(const MoleculeConstraint &other);
    
    virtual ~MoleculeConstraint();
    
    static const char* typeName();
    
    virtual MoleculeConstraint* clone() const=0;

protected:
    MoleculeConstraint& operator=(const MoleculeConstraint &other);
    
    bool operator==(const MoleculeConstraint &other) const;
    bool operator!=(const MoleculeConstraint &other) const;
};

}

SIRE_EXPOSE_CLASS( SireSystem::MoleculeConstraint )

SIRE_END_HEADER

#endif
