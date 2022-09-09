/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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


#include "internalgroupff.h"

#include "SireMaths/line.h"
#include "SireMaths/triangle.h"
#include "SireMaths/torsion.h"

#include "SireBase/property.h"
#include "SireBase/stringproperty.h"
#include "SireBase/propertylist.h"

#include "SireMol/mover.hpp"
#include "SireMol/mgidx.h"

#include "SireFF/detail/atomiccoords3d.h"

#include "SireFF/errors.h"
#include "SireBase/errors.h"
#include "SireError/errors.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/errors.h"

#include "tostring.h"

#include <QDebug>

#include <cstdio>

using namespace SireMM;
using namespace SireMM::detail;
using namespace SireFF;
using namespace SireFF::detail;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireUnits;
using namespace SireStream;

static const RegisterMetaType<InternalGroupFF> r_internalff;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                        const InternalGroupFF &internalff)
{
    writeHeader(ds, r_internalff, 1);

    SharedDataStream sds(ds);

    sds << internalff.mols << internalff.mols0 << internalff.mols1
        << internalff.props
        << internalff.cljgroups
        << internalff.cljgroups0 << internalff.cljgroups1
        << internalff.propmaps
        << static_cast<const G2FF&>(internalff);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                        InternalGroupFF &internalff)
{
    VersionID v = readHeader(ds, r_internalff);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> internalff.mols
            >> internalff.mols0
            >> internalff.mols1
            >> internalff.props
            >> internalff.cljgroups
            >> internalff.cljgroups0
            >> internalff.cljgroups1
            >> internalff.propmaps
            >> static_cast<G2FF&>(internalff);

        internalff._pvt_updateName();
        internalff.calc_14_nrgs = internalff.props.property("calculate14")
                                                  .asABoolean();
    }
    else
        throw version_error(v, "1", r_internalff, CODELOC);

    return ds;
}

/** Constructor */
InternalGroupFF::InternalGroupFF()
                : ConcreteProperty<InternalGroupFF,G2FF>(),
                  FF3D(), InternalPotential(true), calc_14_nrgs(false)
{
    props.setProperty("calculate14", wrap(calc_14_nrgs));
    props.setProperty("combiningRules", wrap("arithmetic"));
}

/** Construct a named internal forcefield */
InternalGroupFF::InternalGroupFF(const QString &name)
                : ConcreteProperty<InternalGroupFF,G2FF>(),
                  FF3D(), InternalPotential(true), calc_14_nrgs(false)
{
    G2FF::setName(name);
    props.setProperty("calculate14", wrap(calc_14_nrgs));
    props.setProperty("combiningRules", wrap("arithmetic"));
}

/** Copy constructor */
InternalGroupFF::InternalGroupFF(const InternalGroupFF &other)
                : ConcreteProperty<InternalGroupFF,G2FF>(other),
                  FF3D(other), InternalPotential(other),
                  mols(other.mols), mols0(other.mols0), mols1(other.mols1),
                  cljgroups(other.cljgroups),
                  cljgroups0(other.cljgroups0),
                  cljgroups1(other.cljgroups1),
                  propmaps(other.propmaps),
                  ffcomponents(other.ffcomponents),
                  props(other.props),
                  calc_14_nrgs(other.calc_14_nrgs)
{}

/** Destructor */
InternalGroupFF::~InternalGroupFF()
{}

/** Copy assignment operator */
InternalGroupFF& InternalGroupFF::operator=(const InternalGroupFF &other)
{
    if (this != &other)
    {
        G2FF::operator=(other);
        FF3D::operator=(other);
        InternalPotential::operator=(other);

        mols = other.mols;
        mols0 = other.mols0;
        mols1 = other.mols1;

        ffcomponents = other.ffcomponents;
        props = other.props;

        cljgroups = other.cljgroups;
        cljgroups0 = other.cljgroups0;
        cljgroups1 = other.cljgroups1;
        propmaps = other.propmaps;
        calc_14_nrgs = other.calc_14_nrgs;
    }

    return *this;
}

/** Comparison operator */
bool InternalGroupFF::operator==(const InternalGroupFF &other) const
{
    return G2FF::operator==(other) and cljgroups == other.cljgroups and
           propmaps == other.propmaps and calc_14_nrgs == other.calc_14_nrgs;
}

/** Comparison operator */
bool InternalGroupFF::operator!=(const InternalGroupFF &other) const
{
    return not this->operator==(other);
}

/** Function used to perform the work of changing the name of this
    forcefield - this renames the component symbols and the molecule group */
void InternalGroupFF::_pvt_updateName()
{
    ffcomponents = Components( this->name() );
    G2FF::_pvt_updateName();
}

/** Turn on the calculate of 1-4 nonbonded terms */
void InternalGroupFF::enable14Calculation()
{
    if (not calc_14_nrgs)
    {
        cljgroups.clear();
        cljgroups0.clear();
        cljgroups1.clear();

        calc_14_nrgs = true;
        CLJFunction::COMBINING_RULES combining_rules = this->combiningRules();

        const ChunkedVector<InternalPotential::Molecule> &mols_array = mols.moleculesByIndex();

        for (int i=0; i<mols.count(); ++i)
        {
            cljgroups.insert( mols_array[i].number(),
                              CLJ14Group(mols_array[i].molecule(), combining_rules, true,
                                         propmaps.value(mols_array[i].number(),
                                                        PropertyMap())) );
        }

        const ChunkedVector<InternalPotential::Molecule> &mols0_array = mols0.moleculesByIndex();

        for (int i=0; i<mols0.count(); ++i)
        {
            cljgroups0.insert( mols0_array[i].number(),
                              CLJ14Group(mols0_array[i].molecule(), combining_rules, true,
                                         propmaps.value(mols0_array[i].number(),
                                                        PropertyMap())) );
        }

        const ChunkedVector<InternalPotential::Molecule> &mols1_array = mols1.moleculesByIndex();

        for (int i=0; i<mols1.count(); ++i)
        {
            cljgroups1.insert( mols1_array[i].number(),
                              CLJ14Group(mols1_array[i].molecule(), combining_rules, true,
                                         propmaps.value(mols1_array[i].number(),
                                                        PropertyMap())) );
        }

        props.setProperty("calculate14", wrap(true));
    }
}

/** Disable calculation of the 1-4 nonbonded terms */
void InternalGroupFF::disable14Calculation()
{
    if (calc_14_nrgs)
    {
        calc_14_nrgs = false;
        cljgroups.clear();

        props.setProperty("calculate14", wrap(false));
    }
}

/** Turn on or off use of arithmetic combining rules when calculating
    the 1-4 nonbonded energy */
void InternalGroupFF::setArithmeticCombiningRules(bool on)
{
    if (on)
        this->setCombiningRules( CLJFunction::ARITHMETIC );
    else
        this->setCombiningRules( CLJFunction::GEOMETRIC );
}

/** Turn on or off use of geometric combining rules when calculating
    the 1-4 nonbonded energy */
void InternalGroupFF::setGeometricCombiningRules(bool on)
{
    if (on)
        this->setCombiningRules( CLJFunction::GEOMETRIC );
    else
        this->setCombiningRules( CLJFunction::ARITHMETIC );
}

/** Return the type of combining rules used when calculating the 1-4
    nonbonded energy */
CLJFunction::COMBINING_RULES InternalGroupFF::combiningRules() const
{
    QString rules = props.property("combiningRules").asAString();

    if (rules == QLatin1String("arithmetic"))
        return CLJFunction::ARITHMETIC;
    else if (rules == QLatin1String("geometric"))
        return CLJFunction::GEOMETRIC;
    else
    {
        qDebug() << "WARNING: COULD NOT INTERPRET COMBINING RULES FROM STRING"
                 << rules << ". USING ARITHMETIC COMBINING RULES";
        return CLJFunction::ARITHMETIC;
    }
}

/** Set the combining rules used when calculating the 1-4 nonbonded energy,
    returning whether or not this changes the forcefield */
bool InternalGroupFF::setCombiningRules(CLJFunction::COMBINING_RULES rules)
{
    CLJFunction::COMBINING_RULES old_rules = this->combiningRules();

    if (rules != old_rules)
    {
        for (QHash<MolNum,CLJ14Group>::iterator it = cljgroups.begin();
             it != cljgroups.end();
             ++it)
        {
            it.value().setCombiningRules(rules);
        }

        switch(rules)
        {
        case CLJFunction::ARITHMETIC:
            props.setProperty("combiningRules", wrap("arithmetic"));
            break;
        case CLJFunction::GEOMETRIC:
            props.setProperty("combiningRules", wrap("geometric"));
            break;
        }

        return false;
    }
    else
        return false;
}

/** Return whether or not arithmetic combining rules are used for the 1-4
    nonbonded energy calculation */
bool InternalGroupFF::usingArithmeticCombiningRules() const
{
    return combiningRules() == CLJFunction::ARITHMETIC;
}

/** Return whether or not geometric combining rules are used for the 1-4
    nonbonded energy calculation */
bool InternalGroupFF::usingGeometricCombiningRules() const
{
    return combiningRules() == CLJFunction::GEOMETRIC;
}

/** Turn on or off the calculation of the 1-4 terms, returning whether
    or not this changes the forcefield */
bool InternalGroupFF::setUse14Calculation(bool on)
{
    if (calc_14_nrgs != on)
    {
        if (on)
            this->enable14Calculation();
        else
            this->disable14Calculation();

        return true;
    }
    else
        return false;
}

/** Return whether or not this forcefield also calculates the
    1-4 nonbonded terms */
bool InternalGroupFF::uses14Calculation() const
{
    return calc_14_nrgs;
}

/** Set the property 'name' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
bool InternalGroupFF::setProperty(const QString &name, const Property &value)
{
    if (name == QLatin1String("combiningRules"))
    {
        QString rules = value.asAString();

        CLJFunction::COMBINING_RULES new_rules;

        if (rules == QLatin1String("arithmetic"))
            new_rules = CLJFunction::ARITHMETIC;
        else if (rules == QLatin1String("geometric"))
            new_rules = CLJFunction::GEOMETRIC;
        else
            throw SireError::invalid_arg( QObject::tr(
                    "Cannot recognise combining rules type from string '%1'. "
                    "Available rules are 'arithmetic' and 'geometric'.")
                        .arg(rules), CODELOC );

        return this->setCombiningRules(new_rules);
    }
    else if (name == QLatin1String("calculate14"))
    {
        return this->setUse14Calculation( value.asABoolean() );
    }
    else
        throw SireBase::missing_property( QObject::tr(
            "InternalGroupFF does not have a property called \"%1\" that "
            "can be changed. Available properties are [ strict ].")
                .arg(name), CODELOC );
}

/** Return the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& InternalGroupFF::property(const QString &name) const
{
    return props.property(name);
}

/** Return whether this forcefield contains the property called 'name' */
bool InternalGroupFF::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return the values of all of the properties of this forcefield */
const Properties& InternalGroupFF::properties() const
{
    return props;
}

/** Calculate the energies in molecules in the passed energy table
    caused by this potential, and add them onto the energies already
    in the energy table (optionally scaled by 'scale_energy') */
void InternalGroupFF::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "InternalGroupFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the energies of molecules in the passed energies table
    caused by the component of this potential represented by
    'symbol', and add them onto the energies already
    in the energy table (optionally scaled by 'scale_energy') */
void InternalGroupFF::energy(EnergyTable &energytable, const Symbol &symbol,
			double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "InternalGroupFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the forces acting on molecules in the passed force table
    caused by this potential, and add them onto the forces already
    in the force table (optionally scaled by 'scale_force') */
void InternalGroupFF::force(ForceTable &forcetable, double scale_force)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the force of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

/** Calculate the forces acting on molecules in the passed force table
    caused by the component of this potential represented by
    'symbol', and add them onto the forces already
    in the force table (optionally scaled by 'scale_force') */
void InternalGroupFF::force(ForceTable &forcetable, const Symbol &symbol,
                            double scale_force)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the force of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

/** Set it that the forcefield must now be recalculate from scratch */
void InternalGroupFF::mustNowRecalculateFromScratch()
{
    //record that the forcefield is dirty
    G2FF::setDirty();
}

/** Internal function used to get a handle on the forcefield components */
const FFComponent& InternalGroupFF::_pvt_components() const
{
    return ffcomponents;
}

/** Recalculate the energy of the current state of this forcefield. This
    will recalculate the energy using the quickest possible route, e.g.
    if will only recalculate the energies of molecules that have changed
    since the last evaluation */
void InternalGroupFF::recalculateEnergy()
{
    // we always recalculate from scratch, because it is much easier!

    //nothing appears to have changed, so lets recalculate
    //everything from scratch
    const int nmols = mols.count();

    const ChunkedVector<InternalPotential::Molecule> &mols_array
                                        = mols.moleculesByIndex();
    const ChunkedVector<InternalPotential::Molecule> &mols0_array
                                        = mols0.moleculesByIndex();
    const ChunkedVector<InternalPotential::Molecule> &mols1_array
                                        = mols1.moleculesByIndex();

    const auto &map0 = mols0.indexesByMolNum();
    const auto &map1 = mols1.indexesByMolNum();

    Energy total_nrg;
    Energy total_nrg0;
    Energy total_nrg1;

    double cnrg = 0;
    double ljnrg = 0;
    double cnrg0 = 0;
    double ljnrg0 = 0;
    double cnrg1 = 0;
    double ljnrg1 = 0;

    for (int i=0; i<nmols; ++i)
    {
        const auto molnum = mols_array[i].number();

        if (map0.contains(molnum) and map1.contains(molnum))
        {
            InternalPotential::calculateEnergy( mols_array[i], total_nrg );
            InternalPotential::calculateEnergy( mols0_array[map0[molnum]], total_nrg0 );
            InternalPotential::calculateEnergy( mols1_array[map1[molnum]], total_nrg1 );

            if (calc_14_nrgs)
            {
                auto it = cljgroups.find(molnum);
                it.value().mustNowRecalculateFromScratch();
                boost::tuple<double,double> nrg = it.value().energy();
                cnrg += nrg.get<0>();
                ljnrg += nrg.get<1>();

                it = cljgroups0.find(molnum);
                it.value().mustNowRecalculateFromScratch();
                nrg = it.value().energy();
                cnrg0 += nrg.get<0>();
                ljnrg0 += nrg.get<1>();

                it = cljgroups1.find(molnum);
                it.value().mustNowRecalculateFromScratch();
                nrg = it.value().energy();
                cnrg1 += nrg.get<0>();
                ljnrg1 += nrg.get<1>();
            }
        }
    }

    total_nrg += Intra14Energy(cnrg, ljnrg);
    total_nrg0 += Intra14Energy(cnrg0, ljnrg0);
    total_nrg1 += Intra14Energy(cnrg1, ljnrg1);

    // we calculate the energy between the two parts as the difference
    // between the total energy of the combined parts minus the
    // energy of each individual part. We use "strict" mode all
    // the time, so the difference must be the inter-group energy
    this->components().setEnergy(*this, total_nrg - total_nrg0 - total_nrg1);

    this->setClean();
}

/** Record that the molecule in 'mol' has been added, using the
    passed property map to get the required forcefield parameters

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void InternalGroupFF::_pvt_added(quint32 group_id,
                                 const PartialMolecule &molecule, const PropertyMap &map)
{
    const auto molnum = molecule.data().number();

    // make sure this molecule is in both groups
    const auto &group0 = this->select(MGIdx(0));
    const auto &group1 = this->select(MGIdx(1));

    if (not (group0.contains(molnum) and group1.contains(molnum)))
    {
        return;
    }

    const auto mol0 = group0[molnum];
    const auto mol1 = group1[molnum];

    if (mol0.selection().intersects(mol1.selection()))
    {
        throw SireMol::duplicate_atom(QObject::tr(
            "One or more atoms of %1:%2 has been added to both groups. You "
            "must be careful to sure that an atom is not added to both groups.")
                .arg(molecule.data().name().value()).arg(molnum.value()),
                    CODELOC);
    }

    mols.add(mol0, map, *this, false);
    mols.add(mol1, map, *this, false);

    mols0.add(mol0, map, *this, false);
    mols1.add(mol1, map, *this, false);

    if (not map.isDefault())
    {
        propmaps.insert(molecule.number(), map);
    }

    if (calc_14_nrgs)
    {
        if (cljgroups.contains(molecule.number()))
        {
            cljgroups[molecule.number()].add(mol0);
            cljgroups[molecule.number()].add(mol1);
        }
        else
        {
            cljgroups.insert(mol0.number(),
                             CLJ14Group(mol0, combiningRules(), true, map));
            cljgroups[molecule.number()].add(mol1);
        }

        if (cljgroups0.contains(molecule.number()))
        {
            cljgroups0[molecule.number()].add(mol0);
        }
        else
        {
            cljgroups0.insert(molecule.number(),
                              CLJ14Group(mol0, combiningRules(), true, map));
        }

        if (cljgroups1.contains(molecule.number()))
        {
            cljgroups1[molecule.number()].add(mol1);
        }
        else
        {
            cljgroups1.insert(molecule.number(),
                              CLJ14Group(mol1, combiningRules(), true, map));
        }
    }
}

/** Record the fact that the molecule 'mol' has been removed from this forcefield */
void InternalGroupFF::_pvt_removed(quint32 group_id,
                                   const PartialMolecule &molecule)
{
    mols.remove(molecule, *this, false);

    if (group_id == 0)
        mols0.remove(molecule, *this, false);
    else
        mols1.remove(molecule, *this, false);

    if (calc_14_nrgs)
    {
        if (cljgroups.contains(molecule.number()))
        {
            cljgroups[molecule.number()].remove(molecule);

            if (cljgroups[molecule.number()].isNull())
            {
                cljgroups.remove(molecule.number());
            }
        }

        if (group_id == 0)
        {
            if (cljgroups0.contains(molecule.number()))
            {
                cljgroups0[molecule.number()].remove(molecule);

                if (cljgroups0[molecule.number()].isNull())
                {
                    cljgroups0.remove(molecule.number());
                }
            }
        }
        else
        {
            if (cljgroups1.contains(molecule.number()))
            {
                cljgroups1[molecule.number()].remove(molecule);

                if (cljgroups1[molecule.number()].isNull())
                {
                    cljgroups1.remove(molecule.number());
                }
            }
        }
    }

    if (not this->contains(molecule.number()))
    {
        propmaps.remove(molecule.number());
    }
}

/** Record that fact that the molecule 'molecule' has been changed

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void InternalGroupFF::_pvt_changed(quint32 group_id,
                                   const SireMol::Molecule &molecule, bool auto_update)
{
    mols.change(molecule, *this, false);
    mols0.change(molecule, *this, false);
    mols1.change(molecule, *this, false);

    if (calc_14_nrgs)
    {
        if (cljgroups.contains(molecule.number()))
        {
            cljgroups[molecule.number()].update(molecule);
        }

        if (cljgroups0.contains(molecule.number()))
        {
            cljgroups0[molecule.number()].update(molecule);
        }

        if (cljgroups1.contains(molecule.number()))
        {
            cljgroups1[molecule.number()].update(molecule);
        }
    }
}

/** Record that the provided list of molecules have changed

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void InternalGroupFF::_pvt_changed(quint32 group_id,
                                   const QList<SireMol::Molecule> &molecules, bool auto_update)
{
    InternalGroupFF::Molecules old_mols = mols;
    QHash<MolNum,CLJ14Group> old_cljgroups = cljgroups;

    InternalGroupFF::Molecules old_mols0 = mols0;
    QHash<MolNum,CLJ14Group> old_cljgroups0 = cljgroups0;

    InternalGroupFF::Molecules old_mols1 = mols1;
    QHash<MolNum,CLJ14Group> old_cljgroups1 = cljgroups1;

    try
    {
        foreach (const SireMol::Molecule &molecule, molecules)
        {
            mols.change(molecule, *this, false);
            mols0.change(molecule, *this, false);
            mols1.change(molecule, *this, false);
        }

        if (calc_14_nrgs)
        {
            foreach (const SireMol::Molecule &molecule, molecules)
            {
                if (cljgroups.contains(molecule.number()))
                {
                    cljgroups[molecule.number()].update(molecule);
                }
                if (cljgroups0.contains(molecule.number()))
                {
                    cljgroups0[molecule.number()].update(molecule);
                }
                if (cljgroups1.contains(molecule.number()))
                {
                    cljgroups1[molecule.number()].update(molecule);
                }
            }
        }
    }
    catch(...)
    {
        mols = old_mols;
        mols0 = old_mols0;
        mols1 = old_mols1;
        cljgroups = old_cljgroups;
        cljgroups0 = old_cljgroups0;
        cljgroups1 = old_cljgroups1;
        throw;
    }
}

/** Record that all of the molecules have been removed */
void InternalGroupFF::_pvt_removedAll(quint32 group_id)
{
    if (group_id == 0)
    {
        mols0.clear();
        cljgroups0.clear();
        mols = mols1;
        cljgroups = cljgroups1;
    }
    else
    {
        mols1.clear();
        cljgroups1.clear();
        mols = mols0;
        cljgroups = cljgroups0;
    }

    this->mustNowRecalculateFromScratch();
}

/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */
bool InternalGroupFF::_pvt_wouldChangeProperties(quint32 group_id, MolNum molnum,
                                                 const PropertyMap &map) const
{
    if (mols.wouldChangeProperties(molnum, map))
        return true;

    else if (cljgroups.contains(molnum))
        return cljgroups[molnum].wouldChangeProperties(map);

    else
        return false;
}

const char* InternalGroupFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalGroupFF>() );
}

InternalGroupFF* InternalGroupFF::clone() const
{
    return new InternalGroupFF(*this);
}

void InternalGroupFF::field(FieldTable &fieldtable, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::field(FieldTable &fieldtable, const Symbol &component,
                       double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::potential(PotentialTable &potentialtable, double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::potential(PotentialTable &potentialtable, const Symbol &component,
                          double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::field(FieldTable &fieldtable, const Probe &probe, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::field(FieldTable &fieldtable, const Symbol &component,
                       const Probe &probe, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::potential(PotentialTable &potentialtable, const Probe &probe,
                           double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

void InternalGroupFF::potential(PotentialTable &potentialtable, const Symbol &component,
                           const Probe &probe, double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalGroupFF has yet to "
                "be implemented."), CODELOC );
}

