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

#include "chargeperturbation.h"
#include "atomcharges.h"
#include "molecule.h"
#include "moleditor.h"
#include "mover.hpp"
#include "selector.hpp"

#include "SireCAS/values.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<ChargePerturbation> r_chgpert;

QDataStream &operator<<(QDataStream &ds,
                                       const ChargePerturbation &chgpert)
{
    writeHeader(ds, r_chgpert, 1);

    ds << static_cast<const Perturbation&>(chgpert);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                       ChargePerturbation &chgpert)
{
    VersionID v = readHeader(ds, r_chgpert);

    if (v == 1)
    {
        ds >> static_cast<Perturbation&>(chgpert);
    }
    else
        throw version_error(v, "1", r_chgpert, CODELOC);

    return ds;
}

/** Constructor - this creates a charge perturbation that
    perturbs from charges in "initial_charge" to charges in
    "final_charge", placing the current charges in "charge",
    and using Perturbation::defaultEquation() to map the
    charges */
ChargePerturbation::ChargePerturbation()
                   : ConcreteProperty<ChargePerturbation,Perturbation>()
{}

/** Construct, using the passed map to find the properties used
    by this perturbation */
ChargePerturbation::ChargePerturbation(const PropertyMap &map)
                   : ConcreteProperty<ChargePerturbation,Perturbation>(map)
{}

/** Construct, using the passed map to find the properties used
    by this perturbation and the passed mapping function to map
    the charges between the states */
ChargePerturbation::ChargePerturbation(const Expression &mapping_function,
                                       const PropertyMap &map)
                   : ConcreteProperty<ChargePerturbation,Perturbation>(mapping_function,
                                                                       map)
{}

/** Copy constructor */
ChargePerturbation::ChargePerturbation(const ChargePerturbation &other)
                   : ConcreteProperty<ChargePerturbation,Perturbation>(other)
{}

/** Destructor */
ChargePerturbation::~ChargePerturbation()
{}

const char* ChargePerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChargePerturbation>() );
}

QString ChargePerturbation::toString() const
{
    return QObject::tr("ChargePerturbation( from %1 to %2 )")
                .arg(propertyMap()["initial_charge"].source(),
                     propertyMap()["final_charge"].source());
}

/** Copy assignment operator */
ChargePerturbation& ChargePerturbation::operator=(const ChargePerturbation &other)
{
    Perturbation::operator=(other);
    return *this;
}

/** Comparison operator */
bool ChargePerturbation::operator==(const ChargePerturbation &other) const
{
    return Perturbation::operator==(other);
}

/** Comparison operator */
bool ChargePerturbation::operator!=(const ChargePerturbation &other) const
{
    return not ChargePerturbation::operator==(other);
}

/** Return the properties required or changed by this perturbation */
QSet<QString> ChargePerturbation::requiredProperties() const
{
    QSet<QString> props;

    PropertyName prop = propertyMap()["charge"];

    if (prop.hasSource())
        props.insert( prop.source() );

    prop = propertyMap()["initial_charge"];

    if (prop.hasSource())
        props.insert( prop.source() );

    prop = propertyMap()["final_charge"];

    if (prop.hasSource())
        props.insert( prop.source() );

    return props;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool ChargePerturbation::wouldChange(const Molecule &molecule,
                                     const Values &values) const
{
    try
    {
        const AtomCharges &initial_chgs = molecule.property(
                                            propertyMap()["initial_charge"] )
                                                .asA<AtomCharges>();

        const AtomCharges &final_chgs = molecule.property(
                                            propertyMap()["final_charge"] )
                                                .asA<AtomCharges>();

        const AtomCharges &chgs = molecule.property(
                                            propertyMap()["charge"] )
                                                .asA<AtomCharges>();

        const Expression &f = this->mappingFunction();
        const Symbol &initial = this->symbols().initial();
        const Symbol &final = this->symbols().final();

        for (CGIdx i(0); i<initial_chgs.nCutGroups(); ++i)
        {
            for (Index j(0); j<initial_chgs.nAtoms(i); ++j)
            {
                CGAtomIdx atomidx(i,j);

                double initial_chg = initial_chgs[atomidx].value();
                double final_chg = final_chgs[atomidx].value();
                double chg = chgs[atomidx].value();

                if (initial_chg != final_chg)
                {
                    Values atom_values = values + (initial == initial_chg) +
                                                  (final == final_chg);

                    if (chg != f(atom_values))
                        return true;
                }
                else if (initial_chg != chg)
                {
                    return true;
                }
            }
        }

        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Perturb the charges in the passed molecule using the reaction
    coordinate(s) in 'values'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ChargePerturbation::perturbMolecule(MolEditor &molecule, const Values &values) const
{
    const AtomCharges &initial_chgs = molecule.property( propertyMap()["initial_charge"] )
                                              .asA<AtomCharges>();

    const AtomCharges &final_chgs = molecule.property( propertyMap()["final_charge"] )
                                            .asA<AtomCharges>();

    AtomCharges chgs(initial_chgs);

    const Expression &f = this->mappingFunction();
    const Symbol &initial = this->symbols().initial();
    const Symbol &final = this->symbols().final();

    for (CGIdx i(0); i<initial_chgs.nCutGroups(); ++i)
    {
        for (Index j(0); j<initial_chgs.nAtoms(i); ++j)
        {
            CGAtomIdx atomidx(i,j);

            double initial_chg = initial_chgs[atomidx].value();
            double final_chg = final_chgs[atomidx].value();

            if (initial_chg != final_chg)
            {
                Values atom_values = values + (initial == initial_chg) +
                                              (final == final_chg);

                chgs.set( atomidx, Charge(f(atom_values)) );
            }
        }
    }

    molecule.setProperty( propertyMap()["charge"].source(), chgs );
}
