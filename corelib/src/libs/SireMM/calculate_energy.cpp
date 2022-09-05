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
#include "SireMM/intraff.h"
#include "SireMM/internalff.h"

#include "SireUnits/units.h"

#include "SireMol/core.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireMol;

using namespace SireBase;

using namespace SireUnits;
using namespace SireUnits::Dimension;

#include <QDebug>

namespace SireMM
{

SIREMM_EXPORT GeneralUnit calculate_energy(ForceFields &ffields)
{
    auto nrgs = ffields.energies();

    double total = ffields.energy().value();

    GeneralUnit nrg;

    for (const auto &ffield : ffields.forceFields())
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
            nrg.addComponent("coulomb_1-4", nrgs[internal.intra14Coulomb()] * kcal_per_mol);
            nrg.addComponent("LJ_1-4", nrgs[internal.intra14LJ()] * kcal_per_mol);
        }
        else if (comps.isA<MultiCLJComponent>())
        {
            const auto &clj = comps.asA<MultiCLJComponent>();

            nrg.addComponent("coulomb_intra", nrgs[clj.coulomb()] * kcal_per_mol);
            nrg.addComponent("LJ_intra", nrgs[clj.lj()] * kcal_per_mol);
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
    ForceFields ffields;

    InternalFF internalff("internal");
    internalff.setStrict(true);
    internalff.add(mol, map);

    IntraFF intraff("intraff");
    intraff.add(mol, map);

    ffields.add(internalff);
    ffields.add(intraff);

    return calculate_energy(ffields);
}

} // end of namespace SireMM
