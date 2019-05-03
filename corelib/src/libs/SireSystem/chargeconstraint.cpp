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

#include "chargeconstraint.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecules.h"
#include "SireMol/molnum.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<ChargeConstraint> r_chgconstraint( MAGIC_ONLY,
                                                    ChargeConstraint::typeName() );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                          const ChargeConstraint &chgconstraint)
{
    writeHeader(ds, r_chgconstraint, 1);
    
    SharedDataStream sds(ds);
    
    sds << chgconstraint.molgroup
        << chgconstraint.prop_map
        << static_cast<const MoleculeConstraint&>(chgconstraint);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          ChargeConstraint &chgconstraint)
{
    VersionID v = readHeader(ds, r_chgconstraint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> chgconstraint.molgroup
            >> chgconstraint.prop_map
            >> static_cast<MoleculeConstraint&>(chgconstraint);
    }
    else
        throw version_error(v, "1", r_chgconstraint, CODELOC);
        
    return ds;
}

/** Constructor */
ChargeConstraint::ChargeConstraint() : MoleculeConstraint()
{}

/** Construct to constrain the charges of the molecules in 'molgroup'
    using the optionally supplied property map to find the properties
    that are needed to calculate and save the charges */
ChargeConstraint::ChargeConstraint(const MoleculeGroup &mgroup, const PropertyMap &map)
                 : MoleculeConstraint(),
                   molgroup(mgroup), prop_map(map)
{}

/** Copy constructor */
ChargeConstraint::ChargeConstraint(const ChargeConstraint &other)
                 : MoleculeConstraint(other),
                   molgroup(other.molgroup), prop_map(other.prop_map)
{}

/** Destructor */
ChargeConstraint::~ChargeConstraint()
{}

const char* ChargeConstraint::typeName()
{
    return "SireSystem::ChargeConstraint";
}

/** Copy assignment operator */
ChargeConstraint& ChargeConstraint::operator=(const ChargeConstraint &other)
{
    if (this != &other)
    {
        molgroup = other.molgroup;
        prop_map = other.prop_map;
        MoleculeConstraint::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool ChargeConstraint::operator==(const ChargeConstraint &other) const
{
    return this == &other or
           (molgroup == other.molgroup and prop_map == other.prop_map and
            MoleculeConstraint::operator==(other));
}

/** Comparison operator */
bool ChargeConstraint::operator!=(const ChargeConstraint &other) const
{
    return not ChargeConstraint::operator==(other);
}

/** Return the molecule group that contains the molecules whose
    charges are being constrained */
const MoleculeGroup& ChargeConstraint::moleculeGroup() const
{
    return molgroup.read();
}

/** Return the property map containing the locations of the properties
    needed to apply this constraint */
const PropertyMap& ChargeConstraint::propertyMap() const
{
    return prop_map;
}

/** Internal function used to match the molecule group help in this
    constraint so that it has the same version of the molecules
    is 'system' */
void ChargeConstraint::updateGroup(const System &system)
{
    if (molgroup.isNull())
        return;
        
    MGNum mgnum = molgroup.read().number();
    
    if (system.contains(mgnum))
        molgroup = system[mgnum];
    else
        molgroup.edit().update(system.molecules());
}
