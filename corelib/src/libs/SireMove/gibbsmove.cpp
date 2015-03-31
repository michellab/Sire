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

/*
#include "gibbsmove.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireStream;
*/
/** Actually perform the Gibbs ensemble move - perform 'nmoves'
    moves, recording simulation statistics if 'record_stats'
    is true */
/*void GibbsMove::move(System &system, int nmoves, bool record_stats)
{
	if (nmoves == 0)
		return;
        
    SaveState old_system_state = SaveState::save(system);

    GibbsMove old_state(*this);
    
    try
    {
        PropertyMap map;
        map.set("coordinates", this->coordinatesProperty());

        for (int i=0; i<nmoves; ++i)
        {
            //save the pre-move system
            System old_system = system;

            //are moving from group0 to group1, or group1 to group0?
            MGIdentifier delete_mgid, add_mgid;
            PropertyName delete_space_name, add_space_name;
            PropertyMap delete_map, add_map;
            Symbol add_nrg, delete_nrg;
        
            if (this->generator().randBool())
            {
                delete_mgid = this->groupID0();
                add_mgid = this->groupID1();
            
                delete_space_name = this->spaceProperty0();
                add_space_name = this->spaceProperty1();
                
                delete_map = map0;
                add_map = map1;
                
                add_nrg = this->
            }
            else
            {
                delete_mgid = this->groupID1();
                add_mgid = this->groupID0();
            
                delete_space_name = this->spaceProperty1();
                add_space_name = this->spaceProperty0();
                
                delete_map = map1;
                add_map = map0;
            }
        
            //get the space in which the old group sits
            const Space &delete_space = system.property(delete_space_name)
                                              .asA<Space>();
                                              
            const Space &add_space = system.property(add_space_name)
                                           .asA<Space>();

            //get the volume of the two simulation boxes
            Volume delete_volume = delete_space.volume();
            Volume add_volume = add_space.volume();
        
            //get the energy of the system before the move
            
        
            //get the group from which a molecule will be deleted
            const MoleculeGroup &delete_group = system[delete_mgid];
        
            //get the group to which the molecule will be added
            const MoleculeGroup &add_group = system[add_mgid];
        
            //pick a molecule at random to delete
            Molecule mol = delete_group[ MolIdx(
                             generator().randInt(delete_group.nMolecules()-1) ) ];

        
            //remove the molecule from the delete_group
            system.remove(mol.number(), delete_mgid);
        
            //now need to randomly translate and rotate this molecule
            Vector rotate_vec = generator().vectorOnSphere();
            Angle rotate_angle = two_pi * generator().rand() * radians;
            
            Vector random_point_in_box = add_space.randomPointInBox( generator() );
            
            Vector mol_center = mol.evaluate().center(add_map);
            
            mol = mol.move().rotate( Quaternion(rotate_vec, rotate_angle),
                                     mol_center, map0 )
                            .translate( random_point_in_box - mol_center )
                            .commit()
        
            //now add this molecule to the new group
            system.add(mol, add_mgid, add_map);
        
            if (collect_stats)
                system.collectStats();
        }
	}
    catch(...)
    {
    	old_system_state.restore(system);
        this->operator=(old_state);
        throw;
    }
}
*/
