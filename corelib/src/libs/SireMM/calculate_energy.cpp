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

#include "calculate_energy.h"

#include "SireMM/interff.h"
#include "SireMM/intergroupff.h"
#include "SireMM/intraff.h"
#include "SireMM/intragroupff.h"
#include "SireMM/internalff.h"
#include "SireMM/intergroupff.h"
#include "SireMM/internalgroupff.h"
#include "SireMM/mmdetail.h"
#include "SireMM/cljshiftfunction.h"

#include "SireMM/intragroupff.h"

#include "SireUnits/units.h"

#include "SireMol/core.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireMol;
using namespace SireVol;

using namespace SireBase;

using namespace SireUnits;
using namespace SireUnits::Dimension;

#include <QDebug>

namespace SireMM
{

MMDetail get_mmdetail (const MoleculeView &mol, const PropertyMap &map)
{
    try
    {
        return mol.data().property(map["forcefield"]).asA<MMDetail>();
    }
    catch(const SireError::exception &e)
    {
        throw SireError::incompatible_error(QObject::tr(
            "Cannot calculate the energy as the molecule %1 does not "
            "contain the necessary metadata describing its forcefield. "
            "You may need to parameterise the molecule first. The "
            "error was %2, %3.")
                .arg(mol.toString())
                .arg(e.what()).arg(e.error()), CODELOC);
    }
}

template<class T>
T get_cljfunc(const MMDetail &mm, SireUnits::Dimension::Length cutoff)
{
    if (not (mm.isAmberStyle() or mm.isOPLS()))
    {
        throw SireError::incompatible_error(QObject::tr(
            "Calculating energies of forcefields that are not Amber- or "
            "OPLS-style is currently not supported. The forcefield style "
            "is %1.").arg(mm.toString()), CODELOC);
    }

    T cljfunc(cutoff);

    if (mm.usesArithmeticCombiningRules())
        cljfunc.setArithmeticCombiningRules(true);
    else if (mm.usesGeometricCombiningRules())
        cljfunc.setGeometricCombiningRules(true);

    return cljfunc;
}

SireUnits::Dimension::Length get_cutoff(const Space &space)
{
    auto cutoff = 15 * angstrom;

    if (space.isPeriodic())
    {
        // work out the maximum size of the box by looking for the periodic
        // image of the cutoff
        auto image = space.getMinimumImage(Vector(0),
                                           Vector(cutoff.value(), 0, 0));

        if (image.x() != cutoff.value())
        {
            // this has been translated
            cutoff = 1.99 * image.x() * angstrom;
            qDebug() << "CUTOFF REDUCED TO" << cutoff.toString();
        }
    }

    return cutoff;
}

SIREMM_EXPORT ForceFields create_forcefield(const MoleculeView &mol,
                                            const PropertyMap &map)
{
    const auto mm = get_mmdetail(mol, map);

    ForceFields ffields;

    InternalFF internalff("internal");
    internalff.setStrict(true);
    internalff.enable14Calculation();

    if (mm.usesArithmeticCombiningRules())
        internalff.setArithmeticCombiningRules(true);
    else if (mm.usesGeometricCombiningRules())
        internalff.setGeometricCombiningRules(true);

    internalff.add(mol, map);

    IntraFF intraff("intraff");
    intraff.setCLJFunction(get_cljfunc<CLJIntraShiftFunction>(mm, 15*angstrom));
    intraff.add(mol, map);

    ffields.add(internalff);
    ffields.add(intraff);

    return ffields;
}

SIREMM_EXPORT ForceFields create_forcefield(const SireMol::Molecules &mols,
                                            const SireBase::PropertyMap &map)
{
    if (mols.isEmpty())
        return ForceFields();

    const auto mm = get_mmdetail(*(mols.constBegin()), map);

    ForceFields ffields;

    // set the space from the first one we can find from the molecules
    const auto space_property = map["space"];
    SpacePtr space = Cartesian();

    for (const auto &mol : mols)
    {
        if (mol.hasProperty(space_property))
        {
            space = mol.data().property(space_property).asA<Space>();
            break;
        }
    }

    InternalFF internalff("internal");
    internalff.setStrict(true);
    internalff.enable14Calculation();

    if (mm.usesArithmeticCombiningRules())
        internalff.setArithmeticCombiningRules(true);
    else if (mm.usesGeometricCombiningRules())
        internalff.setGeometricCombiningRules(true);

    internalff.add(mols, map);

    IntraFF intraff("intraff");
    intraff.setCLJFunction(get_cljfunc<CLJIntraShiftFunction>(mm, 15*angstrom));
    intraff.add(mols, map);

    InterFF interff("interff");
    interff.setCLJFunction(get_cljfunc<CLJShiftFunction>(mm, get_cutoff(*space)));
    interff.setProperty("space", *space);
    interff.add(mols, map);

    ffields.add(internalff);
    ffields.add(intraff);
    ffields.add(interff);

    return ffields;
}

SIREMM_EXPORT ForceFields create_forcefield(const MoleculeView &mol0,
                                            const MoleculeView &mol1,
                                            const PropertyMap &map)
{
    return create_forcefield(Molecules(mol0), Molecules(mol1), map);
}

SIREMM_EXPORT ForceFields create_forcefield(const MoleculeView &mol0,
                                            const Molecules &mols1,
                                            const PropertyMap &map)
{
    return create_forcefield(Molecules(mol0), mols1, map);
}

SIREMM_EXPORT ForceFields create_forcefield(const Molecules &mols0,
                                            const Molecules &mols1,
                                            const PropertyMap &map)
{
    if (mols0.isEmpty() or mols1.isEmpty())
        return ForceFields();

    const auto mm0 = get_mmdetail(*(mols0.constBegin()), map);
    const auto mm1 = get_mmdetail(*(mols1.constBegin()), map);

    if (not mm0.isCompatibleWith(mm1))
    {
        throw SireError::incompatible_error(QObject::tr(
            "Cannot calculate the energy as the molecules have different "
            "types of MM forcefield: %1 versus %2")
                .arg(mm0.toString()).arg(mm1.toString()), CODELOC);
    }

    // set the space from the first one we can find from the molecules
    const auto space_property = map["space"];
    SpacePtr space = Cartesian();
    bool has_property = false;

    for (const auto &mol : mols0)
    {
        if (mol.hasProperty(space_property))
        {
            has_property = true;
            space = mol.data().property(space_property).asA<Space>();
            break;
        }
    }

    if (not has_property)
    {
        for (const auto &mol : mols1)
        {
            if (mol.hasProperty(space_property))
            {
                has_property = true;
                space = mol.data().property(space_property).asA<Space>();
                break;
            }
        }
    }

    ForceFields ffields;

    InterGroupFF interff("interff");
    interff.setCLJFunction(get_cljfunc<CLJShiftFunction>(mm0, get_cutoff(*space)));
    interff.setProperty("space", *space);
    interff.add(mols0, MGIdx(0), map);
    interff.add(mols1, MGIdx(1), map);

    IntraGroupFF intraff("intraff");
    intraff.setCLJFunction(get_cljfunc<CLJIntraShiftFunction>(mm0, 15*angstrom));
    intraff.add(mols0, MGIdx(0), map);
    intraff.add(mols1, MGIdx(1), map);

    InternalGroupFF internalff("internalff");
    internalff.enable14Calculation();

    if (mm0.usesArithmeticCombiningRules())
        internalff.setArithmeticCombiningRules(true);
    else if (mm0.usesGeometricCombiningRules())
        internalff.setGeometricCombiningRules(true);

    internalff.add(mols0, MGIdx(0), map);
    internalff.add(mols1, MGIdx(1), map);

    ffields.add(interff);
    ffields.add(intraff);
    ffields.add(internalff);

    return ffields;
}

SIREMM_EXPORT GeneralUnit calculate_energy(ForceFields &ffields)
{
    auto nrgs = ffields.energies();

    double total = ffields.energy().value();

    GeneralUnit nrg;

    for (const auto &ffield : ffields.forceFields())
    {
        if (ffield->isA<InterFF>() or ffield->isA<InterGroupFF>())
        {
            const auto &clj = ffield.read().components().asA<MultiCLJComponent>();
            nrg.addComponent("coulomb", nrgs[clj.coulomb()] * kcal_per_mol);
            nrg.addComponent("LJ", nrgs[clj.lj()] * kcal_per_mol);
        }
        else if (ffield->isA<IntraFF>() or ffield->isA<IntraGroupFF>())
        {
            const auto &clj = ffield.read().components().asA<MultiCLJComponent>();
            nrg.addComponent("intra_coulomb", nrgs[clj.coulomb()] * kcal_per_mol);
            nrg.addComponent("intra_LJ", nrgs[clj.lj()] * kcal_per_mol);
        }
        else
        {
            const auto &comps = ffield.read().components();

            if (comps.isA<InternalComponent>())
            {
                const auto &internal = comps.asA<InternalComponent>();

                nrg.addComponent("bond", nrgs[internal.bond()] * kcal_per_mol);
                nrg.addComponent("angle", nrgs[internal.angle()] * kcal_per_mol);
                nrg.addComponent("dihedral", nrgs[internal.dihedral()] * kcal_per_mol);
                nrg.addComponent("improper", nrgs[internal.improper()] * kcal_per_mol);
                nrg.addComponent("urey-bradley", nrgs[internal.ureyBradley()] * kcal_per_mol);
                nrg.addComponent("1-4_coulomb", nrgs[internal.intra14Coulomb()] * kcal_per_mol);
                nrg.addComponent("1-4_LJ", nrgs[internal.intra14LJ()] * kcal_per_mol);
            }
            else if (comps.isA<MultiCLJComponent>())
            {
                const auto &clj = comps.asA<MultiCLJComponent>();

                nrg.addComponent("coulomb", nrgs[clj.coulomb()] * kcal_per_mol);
                nrg.addComponent("LJ", nrgs[clj.lj()] * kcal_per_mol);
            }
        }
    }

    double delta = total - nrg.value();

    if (std::abs(delta) > 1e-8)
    {
        nrg.addComponent("other", (total - nrg.value()) * kcal_per_mol);
    }

    return nrg;
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol)
{
    return calculate_energy(mol, PropertyMap());
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol,
                                           const SireBase::PropertyMap &map)
{
    auto ff = create_forcefield(mol, map);
    return calculate_energy(ff);
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols)
{
    return calculate_energy(mols, PropertyMap());
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols,
                                           const SireBase::PropertyMap &map)
{
    auto ff = create_forcefield(mols, map);
    return calculate_energy(ff);
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                           const SireMol::MoleculeView &mol1)
{
    return calculate_energy(mol0, mol1, PropertyMap());
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                           const SireMol::MoleculeView &mol1,
                                           const SireBase::PropertyMap &map)
{
    auto ff = create_forcefield(mol0, mol1, map);
    return calculate_energy(ff);
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                           const SireMol::Molecules &mols1)
{
    return calculate_energy(mol0, mols1, PropertyMap());
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                           const SireMol::Molecules &mols1,
                                           const SireBase::PropertyMap &map)
{
    auto ff = create_forcefield(mol0, mols1, map);
    return calculate_energy(ff);
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols0,
                                           const SireMol::Molecules &mols1)
{
    return calculate_energy(mols0, mols1, PropertyMap());
}

SIREMM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols0,
                                           const SireMol::Molecules &mols1,
                                           const SireBase::PropertyMap &map)
{
    auto ff = create_forcefield(mols0, mols1, map);
    return calculate_energy(ff);
}

} // end of namespace SireMM
