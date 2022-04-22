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
#include "biosimspace.h"

#include "SireBase/getinstalldir.h"

#include "SireError/errors.h"

#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/connectivity.h"
#include "SireMol/mgname.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

using namespace SireBase;
using namespace SireMaths;
using namespace SireMol;
using namespace SireUnits;
using namespace SireVol;

namespace SireIO
{

bool isWater(const Molecule& molecule, const PropertyMap& map)
{
    if (molecule.nAtoms() > 5)
        return false;

    // Get the "element" property from the user map.
    auto elem_prop = map["element"];

    // Zero counters for number of hydrogen and oxygen atoms.
    unsigned num_hydrogen = 0;
    unsigned num_oxygen = 0;

    // Loop over all atoms in the molecule.
    for (int i=0; i<molecule.nAtoms(); ++i)
    {
        const auto atom = molecule.atom(AtomIdx(i));

        try
        {
            const auto element = atom.property<Element>(elem_prop);

            if (element == Element("H"))
                num_hydrogen++;
            else if (element == Element("O"))
                num_oxygen++;
        }
        catch (...)
        {
            // Store the atom name.
            auto name = atom.name().value();

            // Remove all non letter characters.
            name = name.remove(QRegExp("[^a-zA-Z]"));

            // Try to infer the element from the atom name.
            const auto element = Element::biologicalElement(name);

            if (element == Element("H"))
                num_hydrogen++;
            else if (element == Element("O"))
                num_oxygen++;
        }
    }

    if (num_hydrogen == 2 and num_oxygen == 1)
        return true;
    else return false;
}

bool isAmberWater(const Molecule& molecule, const PropertyMap& map)
{
    // Check that this is a water molecule.
    if (not isWater(molecule, map))
        return false;

    // Now check the residue name.
    if (molecule.residue(ResIdx(0)).name().value() != "WAT")
        return false;

    // Now check the atom names.

    // Store the number of atoms.
    const auto num_atoms = molecule.nAtoms();

    // Initialise the atom name template.
    QSet<QString> atom_names;

    // SPC/E or TIP3P.
    if (num_atoms == 3)
    {
        atom_names = QSet<QString>({"O", "H1", "H2"});
    }
    // TIP4P.
    else if (num_atoms == 4)
    {
        atom_names = QSet<QString>({"O", "H1", "H2", "EPW"});
    }
    // TIP5P.
    else if (num_atoms == 5)
    {
        atom_names = QSet<QString>({"O", "H1", "H2", "EP1", "EP2"});
    }

    // Make sure all atom names match the template.
    for (int i=0; i<num_atoms; ++i)
    {
        const auto atom = molecule.atom(AtomIdx(i));

        if (not atom_names.contains(atom.name().value()))
            return false;
    }

    // If we've got this far, then it is an AMBER format water.
    return true;
}

bool isGromacsWater(const Molecule& molecule, const PropertyMap& map)
{
    // Check that this is a water molecule.
    if (not isWater(molecule, map))
        return false;

    // Now check the residue name.
    if (molecule.residue(ResIdx(0)).name().value() != "SOL")
        return false;

    // Now check the atom names.

    // Store the number of atoms.
    const auto num_atoms = molecule.nAtoms();

    // Initialise the atom name template.
    QSet<QString> atom_names;

    // SPC/E or TIP3P.
    if (num_atoms == 3)
    {
        atom_names = QSet<QString>({"OW", "HW1", "HW2"});
    }
    // TIP4P.
    else if (num_atoms == 4)
    {
        atom_names = QSet<QString>({"OW", "HW1", "HW2", "MW"});
    }
    // TIP5P.
    else if (num_atoms == 5)
    {
        atom_names = QSet<QString>({"OW", "HW1", "HW2", "LP1", "LP2"});
    }

    // Make sure all atom names match the template.
    for (int i=0; i<num_atoms; ++i)
    {
        const auto atom = molecule.atom(AtomIdx(i));

        if (not atom_names.contains(atom.name().value()))
            return false;
    }

    // If we've got this far, then it is a GROMACS format water.
    return true;
}

Molecule _pvt_setAmberWater(
    Molecule& molecule,
    const Molecule& water,
    const QString& model,
    bool has_virtual,
    const PropertyMap& map)
{
    // Make the template water molecule editable and renumber it.
    auto edit_mol = water.edit().renumber(molecule.number());

    // Copy across all properties that are unique to the original molecule.
    for (const auto &prop : molecule.propertyKeys())
    {
        if (not molecule.hasProperty(prop))
        {
            edit_mol = edit_mol.setProperty(prop, molecule.property(prop)).molecule();
        }
    }

    // Counter for the number of hydrogen atoms.
    int num_hydrogen = 0;

    // Intialise coordinates vectors for all possible atoms.
    Vector coord_oxygen;
    Vector coord_hydrogen0;
    Vector coord_hydrogen1;

    // Loop over all atoms in the water.
    for (int i=0; i<molecule.nAtoms(); ++i)
    {
        AtomIdx idx(i);
        Element element;

        try
        {
            element = molecule.atom(idx).property<Element>(map["element"]);
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
                coord_hydrogen0 = molecule.atom(idx).property<Vector>(map["coordinates"]);
                num_hydrogen++;
            }
            else
                coord_hydrogen1 = molecule.atom(idx).property<Vector>(map["coordinates"]);
        }

        // Oxygen.
        else if (element == Element("O"))
            coord_oxygen = molecule.atom(idx).property<Vector>(map["coordinates"]);
    }

    // When computing bond potentials for water molecules, AMBER requires
    // that the hydrogen atoms are un-imaged, i.e. close to the oxygen atom.
    // Tools such as gmx trjconv do image water hydrogens, so we need to
    // correct for this.

    // Has the user passed in a "space" property? If so, use this to work
    // get the box dimensions.
    if (map["space"] != "space")
    {
        // Extract the space property.
        auto space = map["space"];

        // A vector to hold the box dimensions.
        Vector box_dims;

        //write out the box dimensions
        if (space.value().isA<PeriodicBox>())
        {
            const auto dims = space.value().asA<PeriodicBox>().dimensions();

            box_dims.setX(dims.x());
            box_dims.setY(dims.y());
            box_dims.setZ(dims.z());
        }
        else if (space.value().isA<TriclinicBox>())
        {
            const auto v0 = space.value().asA<TriclinicBox>().vector0();
            const auto v1 = space.value().asA<TriclinicBox>().vector1();
            const auto v2 = space.value().asA<TriclinicBox>().vector2();

            box_dims.setX(v0.magnitude());
            box_dims.setY(v1.magnitude());
            box_dims.setZ(v2.magnitude());
        }

        // Work out separation vector betwen oxygen and first hydrogen.
        auto sep = coord_hydrogen0 - coord_oxygen;

        // Shift coordinates if separation components exceed half of box.

        // X
        if (std::abs(sep.x()) > 0.5*box_dims.x())
        {
            if (sep.x() < 0)
            {
                coord_hydrogen0.setX(coord_hydrogen0.x() + box_dims.x());
            }
            else
            {
                coord_hydrogen0.setX(coord_hydrogen0.x() - box_dims.x());
            }
        }
        // Y
        if (std::abs(sep.y()) > 0.5*box_dims.y())
        {
            if (sep.y() < 0)
            {
                coord_hydrogen0.setY(coord_hydrogen0.y() + box_dims.y());
            }
            else
            {
                coord_hydrogen0.setY(coord_hydrogen0.y() - box_dims.y());
            }
        }
        // Z
        if (std::abs(sep.z()) > 0.5*box_dims.z())
        {
            if (sep.z() < 0)
            {
                coord_hydrogen0.setZ(coord_hydrogen0.z() + box_dims.z());
            }
            else
            {
                coord_hydrogen0.setZ(coord_hydrogen0.z() - box_dims.z());
            }
        }

        // Now do the same for the second hydrogen.

        sep = coord_hydrogen1 - coord_oxygen;

        // X
        if (std::abs(sep.x()) > 0.5*box_dims.x())
        {
            if (sep.x() < 0)
            {
                coord_hydrogen1.setX(coord_hydrogen1.x() + box_dims.x());
            }
            else
            {
                coord_hydrogen1.setX(coord_hydrogen1.x() - box_dims.x());
            }
        }
        // Y
        if (std::abs(sep.y()) > 0.5*box_dims.y())
        {
            if (sep.y() < 0)
            {
                coord_hydrogen1.setY(coord_hydrogen1.y() + box_dims.y());
            }
            else
            {
                coord_hydrogen1.setY(coord_hydrogen1.y() - box_dims.y());
            }
        }
        // Z
        if (std::abs(sep.z()) > 0.5*box_dims.z())
        {
            if (sep.z() < 0)
            {
                coord_hydrogen1.setZ(coord_hydrogen1.z() + box_dims.z());
            }
            else
            {
                coord_hydrogen1.setZ(coord_hydrogen1.z() - box_dims.z());
            }
        }
    }

    // Replace the atomic coordinates in the template.
    edit_mol = edit_mol.atom(AtomIdx(0)).setProperty(map["coordinates"], coord_oxygen).molecule();
    edit_mol = edit_mol.atom(AtomIdx(1)).setProperty(map["coordinates"], coord_hydrogen0).molecule();
    edit_mol = edit_mol.atom(AtomIdx(2)).setProperty(map["coordinates"], coord_hydrogen1).molecule();

    // Work out coordinates of virtual site(s).
    if (has_virtual)
    {
        // TIP4P
        if (model == "TIP4P")
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

    return edit_mol.commit();
}

Molecule _pvt_setGromacsWater(
    Molecule& molecule,
    const Molecule& water,
    const QString& model,
    bool has_virtual,
    const PropertyMap& map)
{
    // Make the template water molecule editable and renumber it.
    auto edit_mol = water.edit().renumber(molecule.number());

    // Copy across all properties that are unique to the original molecule.
    for (const auto &prop : molecule.propertyKeys())
    {
        if (not molecule.hasProperty(prop))
        {
            edit_mol = edit_mol.setProperty(prop, molecule.property(prop)).molecule();
        }
    }

    // Counter for the number of hydrogen atoms.
    int num_hydrogen = 0;

    // Intialise coordinates vectors for all possible atoms.
    Vector coord_oxygen;
    Vector coord_hydrogen0;
    Vector coord_hydrogen1;

    // Loop over all atoms in the water.
    for (int i=0; i<molecule.nAtoms(); ++i)
    {
        AtomIdx idx(i);
        Element element;

        try
        {
            element = molecule.atom(idx).property<Element>(map["element"]);
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
                coord_hydrogen0 = molecule.atom(idx).property<Vector>(map["coordinates"]);
                num_hydrogen++;
            }
            else
                coord_hydrogen1 = molecule.atom(idx).property<Vector>(map["coordinates"]);
        }

        // Oxygen.
        else if (element == Element("O"))
            coord_oxygen = molecule.atom(idx).property<Vector>(map["coordinates"]);
    }

    // Replace the atomic coordinates in the template.
    edit_mol = edit_mol.atom(AtomIdx(0)).setProperty(map["coordinates"], coord_oxygen).molecule();
    edit_mol = edit_mol.atom(AtomIdx(1)).setProperty(map["coordinates"], coord_hydrogen0).molecule();
    edit_mol = edit_mol.atom(AtomIdx(2)).setProperty(map["coordinates"], coord_hydrogen1).molecule();

    // Work out coordinates of virtual site(s).
    if (has_virtual)
    {
        // TIP4P
        if (model == "TIP4P")
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

    return edit_mol.commit();
}

System setAmberWater(const System& system, const QString& model, const PropertyMap& map)
{
    // Create a new system object.
    System new_system;

    // Create a single molecule group.
    auto molgroup = MoleculeGroup("all");

    // Copy across all system properties.
    for (const auto &prop : system.propertyKeys())
    {
        new_system.setProperty(prop, system.property(prop));
    }

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
    auto water_template = MoleculeParser::read(path + ".prm7", path + ".rst7", map);

    // Extract the water molecule from the template.
    auto template_molecule = water_template[MolIdx(0)].molecule();

    // Loop over all molecules in the system in MolIdx order.
    for (int i=0; i<system.nMolecules(); ++i)
    {
        // Extract the water molecule.
        auto molecule = system.molecule(MolIdx(i)).molecule();

        if (isWater(molecule))
        {
            molecule = _pvt_setAmberWater(molecule,
                template_molecule, _model, has_virtual, map);
        }

        molgroup.add(molecule);
    }

    // Add the molecule group to the system.
    new_system.add(molgroup);

    return new_system;
}

System setGromacsWater(const System& system, const QString& model, const PropertyMap& map)
{
    // Create a new system object.
    System new_system;

    // Create a single molecule group.
    auto molgroup = MoleculeGroup("all");

    // Copy across all system properties.
    for (const auto &prop : system.propertyKeys())
    {
        new_system.setProperty(prop, system.property(prop));
    }

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

    // Loop over all molecules in the system in MolIdx order.
    for (int i=0; i<system.nMolecules(); ++i)
    {
        // Extract the water molecule.
        auto molecule = system.molecule(MolIdx(i)).molecule();

        if (isWater(molecule))
        {
            molecule = _pvt_setGromacsWater(molecule,
                template_molecule, _model, has_virtual, map);
        }

        molgroup.add(molecule);
    }

    // Add the molecule group to the system.
    new_system.add(molgroup);

    return new_system;
}

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
    auto path = models[_model];

    // Load the water model template.
    auto water_template = MoleculeParser::read(path + ".prm7", path + ".rst7", map);

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

        // Apply the new template.
        water = _pvt_setAmberWater(water,
            template_molecule, _model, has_virtual, map);

        // Append the water molecule.
        result.append(water);
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

        // Apply the new template.
        water = _pvt_setGromacsWater(water,
            template_molecule, _model, has_virtual, map);

        // Append the water molecule.
        result.append(water);
    }

    return result;
}

System renumberConstituents(const System& system, unsigned mol_offset)
{
    // Create a new system object.
    System new_system;

    // Create a single molecule group.
    auto molgroup = MoleculeGroup("all");

    // Copy across all system properties.
    for (const auto &prop : system.propertyKeys())
    {
        new_system.setProperty(prop, system.property(prop));
    }

    // Zero the component numbers. These are 1-indexed.
    unsigned num_residues = 1;
    unsigned num_atoms = 1;

    // Work out the cumulative number of constituents up to the
    // molecule offset.
    for (unsigned i=0; i<mol_offset; i++)
    {
        const auto molecule = system.molecule(MolIdx(i)).molecule();

        // Add the molecule.
        molgroup.add(molecule);

        num_residues += molecule.nResidues();
        num_atoms += molecule.nAtoms();
    }

    // Loop over all remaining molecules in the system and renumber
    // their consituents in ascending order.
    for (int i=mol_offset; i<system.nMolecules(); ++i)
    {
        auto molecule = system.molecule(MolIdx(i)).molecule();
        molecule = pvt_renumberConstituents(molecule,
            num_residues, num_atoms);

        molgroup.add(molecule);

        num_residues += molecule.nResidues();
        num_atoms += molecule.nAtoms();
    }

    // Add the molecule group to the system.
    new_system.add(molgroup);

    return new_system;
}

Molecule pvt_renumberConstituents(
    Molecule& molecule,
    const unsigned residue_offset,
    const unsigned atom_offset)
{
    // Make the molecule editable.
    auto edit_mol = molecule.edit();

    // Renumber residues.
    for (int i=0; i<molecule.nResidues(); ++i)
    {
        edit_mol = edit_mol.residue(ResIdx(i)).renumber(ResNum(i+residue_offset)).molecule();
    }

    // Renumber atoms.
    for (int i=0; i<molecule.nAtoms(); ++i)
    {
        edit_mol = edit_mol.atom(AtomIdx(i)).renumber(AtomNum(i+atom_offset)).molecule();
    }

    return edit_mol.commit();
}

System updateAndPreserveOrder(
    const System& system,
    const Molecule& molecule,
    unsigned index)
{
    // Create a new system object.
    System new_system;

    // Create a single molecule group.
    auto molgroup = MoleculeGroup("all");

    // Copy across all system properties.
    for (const auto &prop : system.propertyKeys())
    {
        new_system.setProperty(prop, system.property(prop));
    }

    // Add all of the molecules, preserving the ordering.
    for (unsigned i=0; i<index; ++i)
    {
        auto mol = system.molecule(MolIdx(i)).molecule();
        molgroup.add(mol);
    }

    molgroup.add(molecule);

    for (int i=index+1; i<system.nMolecules(); ++i)
    {
        auto mol = system.molecule(MolIdx(i)).molecule();
        molgroup.add(mol);
    }

    // Add the molecule group to the system.
    new_system.add(molgroup);

    return new_system;
}

System repartitionHydrogenMass(
    const System& system,
    const double factor,
    const unsigned water,
    const PropertyMap& map)
{
    // Create a new system object.
    System new_system;

    // Create a single molecule group.
    auto molgroup = MoleculeGroup("all");

    // Copy across all system properties.
    for (const auto &prop : system.propertyKeys())
    {
        new_system.setProperty(prop, system.property(prop));
    }

    // Loop over all molecules in the system in MolIdx order.
    for (int i=0; i<system.nMolecules(); ++i)
    {
        // Extract the water molecule.
        auto molecule = system.molecule(MolIdx(i)).molecule();

        // Skip water molecules.
        if (water == 0 and isWater(molecule, map))
        {
            molgroup.add(molecule);
            continue;
        }
        // Skip non-water molecules.
        else if (water == 2 and not isWater(molecule, map))
        {
            molgroup.add(molecule);
            continue;
        }

        // This is a perturbable molecule. We need to repartition
        // the mass for both lambda end states.
        if (molecule.hasProperty("is_perturbable"))
        {
            PropertyMap pmap;

            // Lambda = 0 mappings.
            pmap.set("mass", "mass0");
            pmap.set("element", "element0");
            pmap.set("connectivity", "connectivity0");
            pmap.set("coordinates", "coordinates0");

            molecule = repartitionHydrogenMass(
                molecule, factor, water, pmap);

            // Lambda = 1 mappings.
            pmap.set("mass", "mass1");
            pmap.set("element", "element1");
            pmap.set("connectivity", "connectivity1");
            pmap.set("coordinates", "coordinates1");

            molecule = repartitionHydrogenMass(
                molecule, factor, water, pmap);
        }
        else
        {
            molecule = repartitionHydrogenMass(
                molecule, factor, water, map);
        }

        molgroup.add(molecule);
    }

    // Add the molecule group to the system.
    new_system.add(molgroup);

    return new_system;
}

Molecule repartitionHydrogenMass(
    Molecule& molecule,
    const double factor,
    const unsigned water,
    const PropertyMap& map)
{
    // Skip water molecules.
    if (water == 0 and isWater(molecule, map))
    {
        return molecule;
    }
    // Skip non-water molecules.
    else if (water == 2 and not isWater(molecule, map))
    {
        return molecule;
    }

    // Get the name of the element and mass properties.
    const auto elem_prop = map["element"];
    const auto mass_prop = map["mass"];

    if (not molecule.hasProperty(elem_prop))
    {
        throw SireError::incompatible_error(QObject::tr(
            "The molecule doesn't have a \"%1\" property!")
                .arg(elem_prop.toString()), CODELOC);
    }

    if (not molecule.hasProperty(mass_prop))
    {
        throw SireError::incompatible_error(QObject::tr(
            "The molecule doesn't have a \"%1\" property!")
                .arg(mass_prop.toString()), CODELOC);
    }

    // Store the indices of all hydrogen atoms in the system.
    const auto hydrogen = Element("H");
    QList<AtomIdx> hydrogens;
    for (int i=0; i<molecule.nAtoms(); ++i)
    {
        // Extract the atom.
        const auto atom = molecule.atom(AtomIdx(i));

        // Store the index if it is a hydrogen.
        if (atom.property<Element>(elem_prop) == hydrogen)
        {
            hydrogens.append(atom.index());
        }
    }

    // Early exit if this molecule contains no hydrogens.
    if (hydrogens.count() == 0)
        return molecule;

    // Generate the molecular connectivity. (Don't rely on the "connectivity"
    // property.)
    const auto connectivity = Connectivity(molecule, CovalentBondHunter(), map);

    // Compute the initial mass.
    double initial_mass = 0;
    for (const auto &mass : molecule.property(mass_prop).asA<AtomMasses>().toVector())
    {
        initial_mass += mass.value();
    }

    // Make the molecule editable.
    auto edit_mol = molecule.edit();

    // Initialise a list of connections to hydrogen atoms.
    QList<AtomIdx> connections;

    for (const auto &idx : hydrogens)
    {
        // Compute the scaled mass.
        const auto mass = factor * edit_mol.atom(idx).property<Dimension::MolarMass>(mass_prop);
        // Set the new mass.
        edit_mol = edit_mol.atom(idx)
                           .setProperty(mass_prop, mass)
                           .molecule();

        // Store the indices of the atoms that are connected to this hydrogen.
        connections.append(connectivity.connectionsTo(idx).values());
    }

    // Commit the changes.
    molecule = edit_mol.commit();

    // Compute the total adjusted mass.
    double final_mass = 0;
    for (const auto &mass : molecule.property(mass_prop).asA<AtomMasses>().toVector())
    {
        final_mass += mass.value();
    }

    // Work out the delta averaged across the bonded heavy atoms.
    const auto delta_mass = g_per_mol * ((final_mass - initial_mass) / hydrogens.count());

    // Make the molecule editable again.
    edit_mol = molecule.edit();

    // Create a hash mapping the heavy atom index to its current mass.
    QHash<AtomIdx, Dimension::MolarMass> mass_hash;

    // Loop over all connected heavy atoms.
    for (const auto &idx : connections)
    {
        Dimension::MolarMass mass;

        // Use the current mass.
        if (mass_hash.contains(idx))
        {
            mass = mass_hash[idx];
        }
        // Use the initial mass.
        else
        {
            mass =  molecule.atom(idx).property<Dimension::MolarMass>(mass_prop);
        }

        // Reduce the mass.
        mass -= delta_mass;

        // Set the new mass.
        edit_mol = edit_mol.atom(idx)
                           .setProperty(mass_prop, mass)
                           .molecule();

        // Store the updated mass.
        mass_hash[idx] = mass;
    }

    return edit_mol.commit();
}

boost::tuple<System, QHash<MolIdx, MolIdx> > updateCoordinatesAndVelocities(
    const System& system0,
    const System& system1,
    const QHash<MolIdx, MolIdx>& molecule_mapping,
    const bool is_lambda1,
    const PropertyMap& map0,
    const PropertyMap& map1)
{
    // Update coordinates and velocities for the molecules in system0 using
    // the molecules in system1, which may not be in the same order. We assume
    // that both molecules contain unique atom and residue numbers and that
    // atoms retain the same ordering in both systems, i.e. the molecules
    // could be re-ordered, but the atom ordering within those molecules
    // remains the same.

    // Work out the name of the "coordinates" property.
    const auto prop_c0 = map0["coordinates"];
    const auto prop_c1 = map1["coordinates"];

    // Work out the name of the "velocity" property.
    const auto prop_v0 = map0["velocity"];
    const auto prop_v1 = map1["velocity"];

    // Create a local copy of the reference system.
    System new_system = system0;

    // Create a local copy of the molecule mapping.
    QHash<MolIdx, MolIdx> local_mapping(molecule_mapping);

    // Extract the molNums for both systems.
    const auto molNums0 = system0.molNums();
    const auto molNums1 = system1.molNums();

    // A set of molecule indices that have been mapped.
    QSet<int> matched_mols;

    // The mapping is empty, we first need to work out the mapping the
    // molecules in the two systems.
    if (molecule_mapping.isEmpty())
    {
        // Loop over all molecules in system0 in MolIdx order.
        for (int i=0; i<system0.nMolecules(); ++i)
        {
            // Extract the molecule.
            auto molecule0 = system0.molecule(MolIdx(i)).molecule();

            // Get the number of the first atom in the molecule.
            const auto num = molecule0.atom(AtomIdx(0)).number();

            // Loop over all molecules in system1 in MolIdx order.
            for (int j=0; j<system0.nMolecules(); ++j)
            {
                // Make sure this molecule hasn't already been matched.
                if (not matched_mols.contains(j))
                {
                    // Extract the molecule.
                    const auto molecule1 = system1.molecule(MolIdx(j)).molecule();

                    // Extract the first atom from molecule1.
                    const auto atom = molecule1.atom(AtomIdx(0));

                    // The atom numbers match.
                    if (atom.number() == num)
                    {
                        // Extract the molecule.
                        const auto molNum1 = molecule1.number();
                        const auto molecule1 = system1[molNum1].molecule();

                        // Initialise the coordinates property.
                        auto prop_c = prop_c0;

                        // Check whether the molecule is perturbable.
                        if (molecule0.hasProperty("is_perturbable"))
                        {
                            if (is_lambda1) QString
                                prop_c = "coordinates1";
                            else
                                prop_c = "coordinates0";
                        }

                        // Try to update the coordinates property.
                        try
                        {
                            molecule0 = molecule0.edit().setProperty(
                                prop_c, molecule1.property(prop_c1)).molecule().commit();

                            // Update the local mapping.
                            local_mapping[MolIdx(i)] = MolIdx(molNums1.indexOf(molNum1));
                        }
                        catch (...)
                        {
                            throw SireError::incompatible_error(QObject::tr(
                                "Couldn't update coordinates for molecule index '%1'").arg(i), CODELOC);
                        }

                        // Try to update the velocity property. This isn't always present,
                        // so olny try this when the passed system contains the property.
                        if (molecule1.hasProperty(prop_v1))
                        {
                            try
                            {
                                molecule0 = molecule0.edit().setProperty(
                                    prop_v0, molecule1.property(prop_v1)).molecule().commit();
                            }
                            catch (...)
                            {
                                throw SireError::incompatible_error(QObject::tr(
                                    "Couldn't update velocities for molecule index '%1'").arg(i), CODELOC);
                            }
                        }

                        // Update the molecule in the original system.
                        new_system.update(molecule0);

                        // Add the molecule index to the set of matched molecules.
                        matched_mols.insert(j);

                        // We've found a match, so break.
                        break;
                    }
                }
            }
        }

        // Make sure that we've mapped all of the molecules.
        if (matched_mols.count() != system0.nMolecules())
        {
            throw SireError::incompatible_error(QObject::tr(
                "Only matched '%1' molecules, expected '%2'")
                    .arg(matched_mols.count()).arg(system0.nMolecules()), CODELOC);
        }
    }

    // Use the existing molecule mapping.
    else
    {
        // Loop over the mapping.
        for (const auto &molIdx0 : molecule_mapping)
        {
            // Get the molecule index in molecule1.
            const auto molIdx1 = molecule_mapping[molIdx0];

            // Extract the molecules from each system.
            auto molecule0 = system0[molIdx0].molecule();
            const auto molecule1 = system1[molIdx1].molecule();

            // Initialise the coordinates property.
            auto prop_c = prop_c0;

            // Check whether the molecule is perturbable.
            if (molecule0.hasProperty("is_perturbable"))
            {
                if (is_lambda1) QString
                    prop_c = "coordinates1";
                else
                    prop_c = "coordinates0";
            }

            // Try to update the coordinates property.
            try
            {
                molecule0 = molecule0.edit().setProperty(
                    prop_c, molecule1.property(prop_c1)).molecule().commit();
            }
            catch (...)
            {
                throw SireError::incompatible_error(QObject::tr(
                    "Couldn't update coordinates for molecule index '%1'").arg(molIdx0), CODELOC);
            }

            // Try to update the velocity property. This isn't always present,
            // so only try this when the passed system contains the property.
            if (molecule1.hasProperty(prop_v1))
            {
                try
                {
                    molecule0 = molecule0.edit().setProperty(
                        prop_v0, molecule1.property(prop_v1)).molecule().commit();
                }
                catch (...)
                {
                    throw SireError::incompatible_error(QObject::tr(
                        "Couldn't update velocities for molecule index '%1'").arg(molIdx0), CODELOC);
                }
            }

            // Update the molecule in the original system.
            new_system.update(molecule0);
        }
    }

    // Create a tuple containing the updated system and molecule mapping.
    boost::tuple<System, QHash<MolIdx, MolIdx> > retval(new_system, local_mapping);

    return retval;
}

Vector cross(const Vector& v0, const Vector& v1)
{
    double nx = v0.y()*v1.z() - v0.z()*v1.y();
    double ny = v0.z()*v1.x() - v0.x()*v1.z();
    double nz = v0.x()*v1.y() - v0.y()*v1.x();

    return Vector(nx, ny, nz);
}

}
