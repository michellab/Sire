/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIRESYSTEM_CHARGECONSTRAINT_H
#define SIRESYSTEM_CHARGECONSTRAINT_H

#include "moleculeconstraint.h"

#include "SireMol/moleculegroup.h"

#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class ChargeConstraint;
}

QDataStream& operator<<(QDataStream&, const SireSystem::ChargeConstraint&);
QDataStream& operator>>(QDataStream&, SireSystem::ChargeConstraint&);

namespace SireSystem
{

using SireBase::PropertyMap;
using SireBase::PropertyName;

using SireMol::MoleculeGroup;

/** This is the base class of constraints that are used to change 
    the charges on a molecule to match those of an underlying function, 
    e.g. this can be used to recalculate the atomic partial charges of a
    molecule every time it changes conformation, or it can
    be used to modify the charges to correspond to 
    polarisation caused by the environment
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT ChargeConstraint : public MoleculeConstraint
{

friend QDataStream& ::operator<<(QDataStream&, const ChargeConstraint&);
friend QDataStream& ::operator>>(QDataStream&, ChargeConstraint&);

public:
    ChargeConstraint();
    ChargeConstraint(const MoleculeGroup &molgroup,
                     const PropertyMap &map = PropertyMap());
    
    ChargeConstraint(const ChargeConstraint &other);
    
    ~ChargeConstraint();

    static const char* typeName();
    
    const MoleculeGroup& moleculeGroup() const;

    const PropertyMap& propertyMap() const;

protected:    
    ChargeConstraint& operator=(const ChargeConstraint &other);
    
    bool operator==(const ChargeConstraint &other) const;
    bool operator!=(const ChargeConstraint &other) const;

    void updateGroup(const System &system);

	/** The molecule group that contains the molecules
        whose charges are being constrained */
    SireMol::MolGroupPtr molgroup;

    /** The property map used to find the properties that
        are necessary to implement this constraint */
    PropertyMap prop_map;
};

}

SIRE_EXPOSE_CLASS( SireSystem::ChargeConstraint )

SIRE_END_HEADER

#endif
