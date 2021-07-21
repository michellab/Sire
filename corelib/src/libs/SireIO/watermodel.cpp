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

#include "moleculeparser.h"
#include "watermodel.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"

#include "SireSystem/system.h"

using namespace SireBase;
using namespace SireMaths;
using namespace SireMol;

namespace SireIO
{

SelectResult setAmberWater(const SelectResult& molecules, const QString& model, const PropertyMap& map)
{
    // Strip all whitespace from the model name and convert to upper case.
    auto _model = model.simplified().replace(" ", "").toUpper();

    // Create a hash between the allowed model names and their templace files.
    QHash<QString, QString> models;
    models["SPC"]   = getShareDir() + "/templates/water/spce";
    models["SPCE"]  = getShareDir() + "/templates/water/spce";
    models["TIP3P"] = getShareDir() + "/templates/water/tip3p";
    models["TIP4P"] = getShareDir() + "/templates/water/tip4pew";
    models["TIP5P"] = getShareDir() + "/templates/water/tip5p";

    // Make sure the user has passed a valid water model.
    if (not models.contains(_model))
    {
		throw SireError::incompatible_error(QObject::tr(
            "Unsupported AMBER water model '%1'").arg(model), CODELOC);
    }

    // Flag whether the water model has virtual atoms.
    bool has_virtual = false;
    if ((_model == "TIP4P") or (_model == "TIP5P"))
        has_virtual = true;

    // Extract the water model template path.
    _model = models[_model];

    // Load the water model template.
    auto water_template = MoleculeParser::read(_model + ".prm7", _model + ".rst7", map);

    // Extract the water molecule from the template.
    auto template_molecule = water_template[MolIdx(0)].molecule();

    // Make sure that we only operate on water molecules.
    SelectResult waters = molecules.search("water");

    // Create the list of molecules to return.
    QList<ViewsOfMol> result;

    // Loop over all waters in the selection.
    for (const auto &molview : waters.views())
    {
        // Extract the water molecule.
        auto water = molview.molecule();

        // Make the template water molecule editable and renumber it.
        auto edit_mol = template_molecule.edit().renumber();

        // Counter for the number of hydrogen atoms.
        int num_hydrogen = 0;

        // Intialise coordinates vectors for all possible atoms.
        Vector coord_oxygen;
        Vector coord_hydrogen0;
        Vector coord_hydrogen1;

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
                    coord_hydrogen0 = water.atom(idx).property<Vector>(map["coordinates"]);
                    num_hydrogen++;
                }
                else
                    coord_hydrogen1 = water.atom(idx).property<Vector>(map["coordinates"]);
            }

            // Oxygen.
            else if (element == Element("O"))
                coord_oxygen = water.atom(idx).property<Vector>(map["coordinates"]);
        }

        // Replace the atomic coordinates in the template.
        edit_mol = edit_mol.atom(AtomIdx(0)).setProperty(map["coordinates"], coord_oxygen).molecule();
        edit_mol = edit_mol.atom(AtomIdx(1)).setProperty(map["coordinates"], coord_hydrogen0).molecule();
        edit_mol = edit_mol.atom(AtomIdx(2)).setProperty(map["coordinates"], coord_hydrogen1).molecule();

        // Work out coordinates of virtual site(s).
        if (has_virtual)
        {
            // TIP4P
            if (_model == "TIP4P")
            {
                double a = 0.128012065;

                // Expression taken from GROMACS TIP4P topology file.
                // Vsite pos x4 = x1 + a*(x2-x1) + a*(x3-x1)
                // x1 = oxygen, x2 = hydrogen 1, x3 = hydrogn 2

                auto coord_virtual = coord_oxygen + a*(coord_hydrogen0 - coord_oxygen)
                                                  + a*(coord_hydrogen1 - coord_oxygen);

                edit_mol = edit_mol.atom(AtomIdx(3)).setProperty(map["coordinates"], coord_virtual).molecule();
            }

            // TIP5P
            else
            {
                double a = -0.344908262;
                double b = -6.4437903493 / 10.0;

                // Expression taken from GROMACS TIP5P topology file.
                // Vsite pos x4 = x1 + a*(x2-x1) + a*(x3-x1) + b*((x2-x1) x (x3-x1))
                // Vsite pos x5 = x1 + a*(x2-x1) + a*(x3-x1) - b*((x2-x1) x (x3-x1))
                // x1 = oxygen, x2 = hydrogen 1, x3 = hydrogen 2

                auto v0 = coord_hydrogen0 - coord_oxygen;
                auto v1 = coord_hydrogen1 - coord_oxygen;

                auto coord_virtual0 = coord_oxygen + a*(v0 + v1) + b*cross(v0, v1);
                auto coord_virtual1 = coord_oxygen + a*(v0 + v1) - b*cross(v0, v1);

                edit_mol = edit_mol.atom(AtomIdx(3)).setProperty(map["coordinates"], coord_virtual0).molecule();
                edit_mol = edit_mol.atom(AtomIdx(4)).setProperty(map["coordinates"], coord_virtual1).molecule();
            }
        }

        // Append the water molecule.
        result.append(edit_mol.commit());
    }

    return result;
}

SelectResult setGromacsWater(const SelectResult& molecules, const QString& model, const PropertyMap& map)
{
    // Strip all whitespace from the model name and convert to upper case.
    auto _model = model.simplified().replace(" ", "").toUpper();

    // Create a hash between the allowed model names and their templace files.
    QHash<QString, QString> models;
    models["SPC"]   = getShareDir() + "/templates/water/spc";
    models["SPCE"]  = getShareDir() + "/templates/water/spce";
    models["TIP3P"] = getShareDir() + "/templates/water/tip3p";
    models["TIP4P"] = getShareDir() + "/templates/water/tip4p";
    models["TIP5P"] = getShareDir() + "/templates/water/tip5p";

    // Make sure the user has passed a valid water model.
    if (not models.contains(_model))
    {
		throw SireError::incompatible_error(QObject::tr(
            "Unsupported AMBER water model '%1'").arg(model), CODELOC);
    }

    // Flag whether the water model has virtual atoms.
    bool has_virtual = false;
    if ((_model == "TIP4P") or (_model == "TIP5P"))
        has_virtual = true;

    // Extract the water model template path.
    auto path = models[_model];

    // Load the water model template.
    auto water_template = MoleculeParser::read(path + ".grotop", path + ".gro87", map);

    // Extract the water molecule from the template.
    auto template_molecule = water_template[MolIdx(0)].molecule();

    // Make sure that we only operate on water molecules.
    SelectResult waters = molecules.search("water");

    // Creat the list of molecules to return.
    QList<ViewsOfMol> result;

    // Loop over all waters in the selection.
    for (const auto &molview : waters.views())
    {
        // Extract the water molecule.
        auto water = molview.molecule();

        // Make the template water molecule editable and renumber it.
        auto edit_mol = template_molecule.edit().renumber();

        // Counter for the number of hydrogen atoms.
        int num_hydrogen = 0;

        // Intialise coordinates vectors for all possible atoms.
        Vector coord_oxygen;
        Vector coord_hydrogen0;
        Vector coord_hydrogen1;

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
                    coord_hydrogen0 = water.atom(idx).property<Vector>(map["coordinates"]);
                    num_hydrogen++;
                }
                else
                    coord_hydrogen1 = water.atom(idx).property<Vector>(map["coordinates"]);
            }

            // Oxygen.
            else if (element == Element("O"))
                coord_oxygen = water.atom(idx).property<Vector>(map["coordinates"]);
        }

        // Replace the atomic coordinates in the template.
        edit_mol = edit_mol.atom(AtomIdx(0)).setProperty(map["coordinates"], coord_oxygen).molecule();
        edit_mol = edit_mol.atom(AtomIdx(1)).setProperty(map["coordinates"], coord_hydrogen0).molecule();
        edit_mol = edit_mol.atom(AtomIdx(2)).setProperty(map["coordinates"], coord_hydrogen1).molecule();

        // Work out coordinates of virtual site(s).
        if (has_virtual)
        {
            // TIP4P
            if (_model == "TIP4P")
            {
                double a = 0.128012065;

                // Expression taken from GROMACS TIP4P topology file.
                // Vsite pos x4 = x1 + a*(x2-x1) + a*(x3-x1)
                // x1 = oxygen, x2 = hydrogen 1, x3 = hydrogn 2

                auto coord_virtual = coord_oxygen + a*(coord_hydrogen0 - coord_oxygen)
                                                  + a*(coord_hydrogen1 - coord_oxygen);

                edit_mol = edit_mol.atom(AtomIdx(3)).setProperty(map["coordinates"], coord_virtual).molecule();
            }

            // TIP5P
            else
            {
                double a = -0.344908262;
                double b = -6.4437903493 / 10.0;

                // Expression taken from GROMACS TIP5P topology file.
                // Vsite pos x4 = x1 + a*(x2-x1) + a*(x3-x1) + b*((x2-x1) x (x3-x1))
                // Vsite pos x5 = x1 + a*(x2-x1) + a*(x3-x1) - b*((x2-x1) x (x3-x1))
                // x1 = oxygen, x2 = hydrogen 1, x3 = hydrogen 2

                auto v0 = coord_hydrogen0 - coord_oxygen;
                auto v1 = coord_hydrogen1 - coord_oxygen;

                auto coord_virtual0 = coord_oxygen + a*(v0 + v1) + b*cross(v0, v1);
                auto coord_virtual1 = coord_oxygen + a*(v0 + v1) - b*cross(v0, v1);

                edit_mol = edit_mol.atom(AtomIdx(3)).setProperty(map["coordinates"], coord_virtual0).molecule();
                edit_mol = edit_mol.atom(AtomIdx(4)).setProperty(map["coordinates"], coord_virtual1).molecule();
            }
        }

        // Append the water molecule.
        result.append(edit_mol.commit());
    }

    return result;
}

Vector cross(const Vector& v0, const Vector& v1)
{
    double nx = v0.y()*v1.z() - v0.z()*v1.y();
    double ny = v0.z()*v1.x() - v0.x()*v1.z();
    double nz = v0.x()*v1.y() - v0.y()*v1.x();

    return Vector(nx,ny,nz);
}

}
