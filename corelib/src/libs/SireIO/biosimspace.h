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

#ifndef SIREIO_BIOSIMSPACE_H
#define SIREIO_BIOSIMSPACE_H

#include <boost/tuple/tuple.hpp>

#include "sireglobal.h"

#include "SireBase/propertymap.h"

#include "SireMaths/vector.h"

#include "SireMol/select.h"

SIRE_BEGIN_HEADER

using namespace SireBase;
using namespace SireMaths;
using namespace SireMol;
using namespace SireSystem;

namespace SireIO
{
    //! Test whether the passed water molecule matches standard water
    /*! topologies.

        \param molecule
            The molecule to test.

        \param map
            A dictionary of user-defined molecular property names.

        \retval is_water
            Whether the molecule is a water.
     */
    SIREIO_EXPORT bool isWater(
            const Molecule& molecule,
            const PropertyMap& map = PropertyMap());

    //! Test whether the passed water molecule matches standard AMBER
    /*! format water topologies.

        \param molecule
            The molecule to test.

        \param map
            A dictionary of user-defined molecular property names.

        \retval is_water
            Whether the molecule is an AMBER format water.
     */
    SIREIO_EXPORT bool isAmberWater(
            const Molecule& molecule,
            const PropertyMap& map = PropertyMap());

    //! Test whether the passed water molecule matches standard GROMACS
    /*! format water topologies.

        \param molecule
            The molecule to test.

        \param map
            A dictionary of user-defined molecular property names.

        \retval is_water
            Whether the molecule is a GROMACS format water.
     */
    SIREIO_EXPORT bool isGromacsWater(
            const Molecule& molecule,
            const PropertyMap& map = PropertyMap());

    Molecule _pvt_setAmberWater(
            Molecule& molecule,
            const Molecule& water,
            const QString& model,
            bool has_virtual,
            const PropertyMap& map = PropertyMap());

    Molecule _pvt_setGromacsWater(
            Molecule& molecule,
            const Molecule& water,
            const QString& model,
            bool has_virtual,
            const PropertyMap& map = PropertyMap());

    //! Set all water molecules in the passed system to the appropriate AMBER
    /*! format topology.

        \param system
            The molecular system of interest.

        \param model
            The name of the water model.

        \param map
            A dictionary of user-defined molecular property names.

        \retval system
            The system with updated water topology.
     */
    SIREIO_EXPORT System setAmberWater(
            const System& system,
            const QString& model,
            const PropertyMap& map = PropertyMap());

    //! Set all water molecules in the passed system to the appropriate GROMACS
    /*! format topology.

        \param system
            The molecular system of interest.

        \param model
            The name of the water model.

        \param map
            A dictionary of user-defined molecular property names.

        \retval system
            The system with updated water topology.
     */
    SIREIO_EXPORT System setGromacsWater(
            const System& system,
            const QString& model,
            const PropertyMap& map = PropertyMap());

    //! Set all water molecules in the passed system to the appropriate AMBER
    /*! format topology.

        \param system
            The molecular system of interest.

        \param model
            The name of the water model.

        \param map
            A dictionary of user-defined molecular property names.

        \retval system
            The system with updated water topology.
     */
    SIREIO_EXPORT SelectResult setAmberWater(
            const SelectResult& molecules,
            const QString& model,
            const PropertyMap& map = PropertyMap());

    //! Set all water molecules in the passed system to the appropriate GROMACS
    /*! format topology.

        \param system
            The molecular system of interest.

        \param model
            The name of the water model.

        \param map
            A dictionary of user-defined molecular property names.

        \retval system
            The system with updated water topology.

     */
    SIREIO_EXPORT SelectResult setGromacsWater(
            const SelectResult& molecules,
            const QString& model,
            const PropertyMap& map = PropertyMap());

    //! Renumber the constituents of a system (residues and atoms) so that
    /*! they are unique and are in ascending order.

        \param system
            The molecular system of interest.

        \param mol_offset
            The index of the molecule at which to begin renumbering.

        \retval system
            The system with renumbered constituents.
     */
    SIREIO_EXPORT System renumberConstituents(
            const System& system,
            unsigned mol_offset=0);

    Molecule pvt_renumberConstituents(
            Molecule& molecule,
            const unsigned residue_offset,
            const unsigned atom_offset);

    //! Update a molecule in the system with a different UUID while
    /*! preserving the molecular ordering. Normally we would need to
        delete and re-add the molecule, which would place it at the
        end, even if the MolNum was unchanged.

        \param system
            The molecular system of interest.

        \param molecule
            The updated molecule.

        \param index
            The index of the molecule in the system.

        \retval system
            The system with renumbered constituents.
     */
    SIREIO_EXPORT System updateAndPreserveOrder(
            const System& system,
            const Molecule& molecule,
            unsigned index);

    //! Redistribute mass of heavy atoms connected to bonded hydrogens into
    /*! the hydrogen atoms. This allows use of larger simulation integration
        time steps without encountering instabilities related to high-frequency
        hydrogen motion.

        \param system
            The molecular system of interest.

        \param factor
            The repartitioning scale factor. Hydrogen masses are scaled by
            this amount.

        \param water
            Whether to repartiotion masses for water molecules:
            0 = yes, 1 = no, 2 = only water molecules.

        \param map
            A dictionary of user-defined molecular property names.

        \retval system
            The system with repartitioned hydrogen mass.
     */
    SIREIO_EXPORT System repartitionHydrogenMass(
            const System& system,
            const double factor=4,
            const unsigned water=0,
            const PropertyMap& map = PropertyMap());

    //! Redistribute mass of heavy atoms connected to bonded hydrogens into
    /*! the hydrogen atoms. This allows use of larger simulation integration
        time steps without encountering instabilities related to high-frequency
        hydrogen motion.

        \param molecule
            The molecule of interest.

        \param factor
            The repartitioning scale factor. Hydrogen masses are scaled by
            this amount.

        \param water
            Whether to repartiotion masses for water molecules:
            0 = yes, 1 = no, 2 = only water molecules.

        \param map
            A dictionary of user-defined molecular property names.

        \retval system
            The system with repartitioned hydrogen mass.
     */
    SIREIO_EXPORT Molecule repartitionHydrogenMass(
            Molecule& molecule,
            const double factor=4,
            const unsigned water=0,
            const PropertyMap& map = PropertyMap());

    //! Update the coordinates and velocities of system0 with those from
    /*! system1.

        \param system0
            The reference system.

        \param system1
            The updated system, where molecules may not be in the same order.

        \param map0
            A dictionary of user-defined molecular property names for system0.

        \param map1
            A dictionary of user-defined molecular property names for system1.

        \retval system, mapping
            The system with updated coordinates and velocities and a mapping
            between the molecule indices in both systems.
     */
    SIREIO_EXPORT boost::tuple<System, QHash<MolIdx, MolIdx> > updateCoordinatesAndVelocities(
            const System& system0,
            const System& system1,
            const QHash<MolIdx, MolIdx>& molecule_mapping,
            const bool is_lambda1=false,
            const PropertyMap& map0 = PropertyMap(),
            const PropertyMap& map1 = PropertyMap());

    Vector cross(const Vector& v0, const Vector& v1);
}

SIRE_EXPOSE_FUNCTION( SireIO::isWater )
SIRE_EXPOSE_FUNCTION( SireIO::isAmberWater )
SIRE_EXPOSE_FUNCTION( SireIO::isGromacsWater )
SIRE_EXPOSE_FUNCTION( SireIO::setAmberWater )
SIRE_EXPOSE_FUNCTION( SireIO::setGromacsWater )
SIRE_EXPOSE_FUNCTION( SireIO::renumberConstituents )
SIRE_EXPOSE_FUNCTION( SireIO::updateAndPreserveOrder )
SIRE_EXPOSE_FUNCTION( SireIO::repartitionHydrogenMass )
SIRE_EXPOSE_FUNCTION( SireIO::updateCoordinatesAndVelocities )

SIRE_END_HEADER

#endif
