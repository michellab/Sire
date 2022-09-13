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

#include "getcoordsarray.h"

#include "SireMol/core.h"
#include "SireMol/molidx.h"

#include "SireBase/parallel.h"

#include "SireUnits/units.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireMaths;
using namespace SireSystem;

namespace SireIO
{
    QVector<float>
    _getCoordsArray(const Selector<Atom> &atoms,
                    const float scl,
                    const PropertyName &coords_property)
    {
        const int natoms = atoms.count();

        QVector<float> ret(3*natoms, 0.0);
        float *ret_data = ret.data();

        for (int i=0; i<natoms; ++i)
        {
            const int idx = 3 * i;

            const auto &atom = atoms(i);

            const Vector &coords = atom.property<Vector>(coords_property);

            ret_data[idx] = coords.x() * scl;
            ret_data[idx+1] = coords.y() * scl;
            ret_data[idx+2] = coords.z() * scl;
        }

        return ret;
    }


    SIREIO_EXPORT QVector<float>
    getCoordsArray(const SireMol::MoleculeView &mol,
                   const SireUnits::Dimension::Length &to_unit,
                   const SireBase::PropertyMap &map)
    {
        return _getCoordsArray(mol.atoms(),
                               SireUnits::angstrom.to(to_unit),
                               map["coordinates"]);
    }


    SIREIO_EXPORT QVector<float>
    getCoordsArray(const SireMol::MoleculeGroup &mols,
                   const SireUnits::Dimension::Length &to_unit,
                   const SireBase::PropertyMap &map)
    {
        const auto coords_property = map["coordinates"];
        const float scl = SireUnits::angstrom.to(to_unit);

        const int nmols = mols.nMolecules();

        if (nmols == 0)
        {
            return QVector<float>();
        }
        else if (nmols == 1)
        {
            return _getCoordsArray(mols[MolIdx(0)].atoms(),
                                   scl, coords_property);
        }

        QVector< QVector<float> > ret(nmols);
        QVector<float> *ret_data = ret.data();

        if (should_run_in_parallel(nmols, map))
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,nmols),
                               [&](const tbb::blocked_range<int> &r){

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    ret_data[i] = _getCoordsArray(mols[MolIdx(i)].atoms(),
                                                  scl, coords_property);
                }
            });
        }
        else
        {
            for (int i=0; i<nmols; ++i)
            {
                ret_data[i] = _getCoordsArray(mols[MolIdx(i)].atoms(),
                                              scl, coords_property);
            }
        }

        int n = 0;

        for (int i=0; i<nmols; ++i)
        {
            n += ret_data[i].count();
        }

        QVector<float> r(n, 0.0);
        float *r_data = r.data();

        for (int i=0; i<nmols; ++i)
        {
            std::memcpy(r_data, ret_data[i].constData(),
                        ret_data[i].count() * sizeof(float));

            r_data += ret_data[i].count();
        }

        return r;
    }

    SIREIO_EXPORT QVector<float>
    getCoordsArray(const SireSystem::System &system,
                   const SireUnits::Dimension::Length &to_unit,
                   const SireBase::PropertyMap &map)
    {
        MoleculeGroup tmp("tmp");

        for (int i=0; i<system.nMolecules(); ++i)
        {
            tmp.add(system[MolIdx(i)]);
        }

        return getCoordsArray(tmp, to_unit, map);
    }
}
