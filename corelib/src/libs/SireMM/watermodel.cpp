/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Lester Hedges
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

#include "watermodel.h"

#include "SireMol/atomelements.h"
#include "SireMol/moleditor.h"

using namespace SireMol;

namespace SireMM
{

SelectResult SIREMM_EXPORT setAmberWater(const SelectResult& molecules, const PropertyMap& map)
{
    // Make sure that we only operate on water molecules.
    SelectResult waters = molecules.search("water");

    QList<ViewsOfMol> result;

    // Loop over all waters in the selection.
    for (const auto &molview : waters.views())
    {
        // Extract the water molecule and make it editable.
        auto water = molview.molecule().edit();

        // Counter for the number of hydrogen atoms.
        int num_hydrogen = 0;

        // Loop over all atoms in the water.
        for (int i=0; i<water.nAtoms(); ++i)
        {
            AtomIdx idx(i);
            Element element;

            try
            {
                element = water.atom(idx).property<Element>(map["element"]);
            }
            catch (...)
            {
                continue;
            }

            // Hydrogen.
            if (element == Element("H"))
            {
                if (num_hydrogen == 0)
                {
                    water = water.atom(idx).rename(AtomName("H1")).molecule().commit();
                    num_hydrogen++;
                }
                else
                    water = water.atom(idx).rename(AtomName("H2")).molecule().commit();
            }

            // Oxygen.
            else if (element == Element("O"))
                water = water.atom(idx).rename(AtomName("O")).molecule().commit();

            // Dummy.
            else if (element == Element("X"))
                water = water.atom(idx).rename(AtomName("EPW")).molecule().commit();
        }

        // Append the updated water molecule.
        result.append(water.commit());
    }

    return result;
}

SelectResult SIREMM_EXPORT setGromacsWater(const SelectResult& molecules, const PropertyMap& map)
{
    // Make sure that we only operate on water molecules.
    SelectResult waters = molecules.search("water");

    QList<ViewsOfMol> result;

    // Loop over all waters in the selection.
    for (const auto &molview : waters.views())
    {
        // Extract the water molecule and make it editable.
        auto water = molview.molecule().edit();

        // Counter for the number of hydrogen atoms.
        int num_hydrogen = 0;

        // Loop over all atoms in the water.
        for (int i=0; i<water.nAtoms(); ++i)
        {
            AtomIdx idx(i);
            Element element;

            try
            {
                element = water.atom(idx).property<Element>(map["element"]);
            }
            catch (...)
            {
                continue;
            }

            // Hydrogen.
            if (element == Element("H"))
            {
                if (num_hydrogen == 0)
                {
                    water = water.atom(idx).rename(AtomName("HW1")).molecule().commit();
                    num_hydrogen++;
                }
                else
                    water = water.atom(idx).rename(AtomName("HW2")).molecule().commit();
            }

            // Oxygen.
            else if (element == Element("O"))
                water = water.atom(idx).rename(AtomName("OW")).molecule().commit();

            // Dummy.
            else if (element == Element("X"))
                water = water.atom(idx).rename(AtomName("MW")).molecule().commit();
        }

        // Append the updated water molecule.
        result.append(water.commit());
    }

    return result;
}

}
