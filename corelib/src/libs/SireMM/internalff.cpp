/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "internalff.h"

#include "SireMaths/line.h"
#include "SireMaths/triangle.h"
#include "SireMaths/torsion.h"

#include "SireBase/property.h"
#include "SireBase/stringproperty.h"
#include "SireBase/propertylist.h"

#include "SireMol/mover.hpp"

#include "SireFF/detail/atomiccoords3d.h"

#include "SireFF/errors.h"
#include "SireBase/errors.h"
#include "SireError/errors.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

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

////////
//////// Fully instantiate template classes
////////

namespace SireFF
{
    namespace detail
    {
        template class FFMolecule<InternalPotential>;
        template class FFMolecules<InternalPotential>;
        template class ChangedMolecule<InternalPotential::Molecule>;
    }
}

//////// Parameter names

PropertyName BondParameterName::bond_param( "bond", Property::null() );
PropertyName AngleParameterName::angle_param( "angle", Property::null() );
PropertyName DihedralParameterName::dihedral_param( "dihedral", Property::null() );
PropertyName ImproperParameterName::improper_param( "improper", Property::null() );
PropertyName UreyBradleyParameterName::ub_param( "Urey-Bradley", Property::null() );
PropertyName StretchStretchParameterName::ss_param( "stretch-stretch", 
                                                    Property::null() );
PropertyName StretchBendParameterName::sb_param( "stretch-bend", 
                                                 Property::null() );
PropertyName BendBendParameterName::bb_param( "bend-bend", Property::null() );
PropertyName StretchBendTorsionParameterName::sbt_param( "stretch-bend-torsion", 
                                                         Property::null() );

////////
//////// Implementation of InternalPotential
////////

/** Constructor - if 'isstrict' is true then this only calculates
    the energies and forces of parameters that only involve selected
    atoms. If 'isstrict' is false, then parameters that involve
    at least one selected atom are used. */
InternalPotential::InternalPotential(bool strict) : isstrict(strict)
{}

/** Copy constructor */
InternalPotential::InternalPotential(const InternalPotential &other)
                  : isstrict(other.isstrict)
{}

/** Destructor */
InternalPotential::~InternalPotential()
{}

/** Copy assignment operator */
InternalPotential& InternalPotential::operator=(const InternalPotential &other)
{
    isstrict = other.isstrict;
    return *this;
}

InternalPotential::ParameterNames InternalPotential::param_names;
InternalSymbols InternalPotential::internal_symbols;

/** Return all of the symbols used in the internal energy functions */
const InternalSymbols& InternalPotential::symbols()
{
    return internal_symbols;
}

/** Return the names of all of the properties used to store the 
    parameters for this potential */
const InternalPotential::ParameterNames& InternalPotential::parameters()
{
    return param_names;
}

/** Return all of the required parameters for the molecule 'molecule',
    using the supplied map to find the required properties
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InternalPotential::Parameters
InternalPotential::getParameters(const PartialMolecule &molecule,
                                 const PropertyMap &map) const
{
    return InternalPotential::Parameters(molecule,
                                         map[parameters().coordinates()],
                                         map[parameters().bond()],
                                         map[parameters().angle()],
                                         map[parameters().dihedral()],
                                         map[parameters().improper()],
                                         map[parameters().ureyBradley()],
                                         map[parameters().stretchStretch()],
                                         map[parameters().stretchBend()],
                                         map[parameters().bendBend()],
                                         map[parameters().stretchBendTorsion()],
                                         isstrict);
}

static inline bool changed(const PartialMolecule &old_molecule,
                           const PartialMolecule &new_molecule,
                           const PropertyName &property)
{
    return old_molecule.version(property) != new_molecule.version(property);
}

/** Return updated parameters for the molecule that has changed from 
    'old_molecule' to 'new_molecule', with the map 'map' used to find
    the properties of these parameters in both versions of the molecule.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InternalPotential::Parameters
InternalPotential::updateParameters(const InternalPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &map) const
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just return completely
        //new parameters
        return this->getParameters(new_molecule, map);
        
    //we can only update things if just the coordinates have changed
    if ( changed(old_molecule, new_molecule, map[parameters().bond()]) or
         changed(old_molecule, new_molecule, map[parameters().angle()]) or
         changed(old_molecule, new_molecule, map[parameters().dihedral()]) or
         changed(old_molecule, new_molecule, map[parameters().improper()]) or
         changed(old_molecule, new_molecule, map[parameters().ureyBradley()]) or
         changed(old_molecule, new_molecule, map[parameters().stretchStretch()]) or
         changed(old_molecule, new_molecule, map[parameters().stretchBend()]) or
         changed(old_molecule, new_molecule, map[parameters().bendBend()]) or
         changed(old_molecule, new_molecule, map[parameters().stretchBendTorsion()]) )
    {
        //the internal parameters have changed, so we need to update
        //everything
        return this->getParameters(new_molecule, map);
    }
    
    const PropertyName &coords_property = map[parameters().coordinates()];
    
    if ( changed(old_molecule, new_molecule, coords_property) )
    {
        //just the coordinates have changed - just update these
        Parameters new_params = old_params;
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule.molecule(),
                                                        coords_property) );
                                          
        return new_params;
    }
    else
        //nothing important has changed
        return old_params;
}

static inline bool changed(const PartialMolecule &old_molecule,
                           const PartialMolecule &new_molecule,
                           const PropertyName &old_property,
                           const PropertyName &new_property)
{
    return old_property != new_property or
           old_molecule.version(old_property) != 
                                    new_molecule.version(new_property);
}

/** Return updated parameters for the molecule that has changed from 
    'old_molecule' to 'new_molecule', with the map 'old_map' used to find
    the properties of these parameters in the old version of the molecule,
    and 'new_map' used to find the properties in the new version
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InternalPotential::Parameters
InternalPotential::updateParameters(const InternalPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &old_map,
                                    const PropertyMap &new_map) const
{
    if (old_molecule.selection() != new_molecule.selection())
        //we need to rebuild the parameters from scratch
        return this->getParameters(new_molecule, new_map);
        
    if ( changed(old_molecule, new_molecule,
                 old_map[parameters().bond()], new_map[parameters().bond()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().angle()], new_map[parameters().angle()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().dihedral()], new_map[parameters().dihedral()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().improper()], new_map[parameters().improper()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().ureyBradley()], 
                 new_map[parameters().ureyBradley()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().stretchStretch()], 
                 new_map[parameters().stretchStretch()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().stretchBend()], 
                 new_map[parameters().stretchBend()]) or
         changed(old_molecule, new_molecule,
                 old_map[parameters().stretchBendTorsion()], 
                 new_map[parameters().stretchBendTorsion()]) )
    {
        //some of the parameters have changed - we have to rebuild from scratch
        return this->getParameters(new_molecule, new_map);
    }

    if ( changed(old_molecule, new_molecule,
                 old_map[parameters().coordinates()],
                 new_map[parameters().coordinates()]) )
    {
        //only the coordinates have changed
        Parameters new_params = old_params;
        
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule.molecule(),
                                              new_map[parameters().coordinates()]) );
                                              
        return new_params;
    }
    else
        //nothing relavant has changed
        return old_params;
}

/** Return a fully parameterised version 'molecule', using the supplied
    property map to find the properties containing the required 
    parameters.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InternalPotential::Molecule
InternalPotential::parameterise(const PartialMolecule &molecule,
                                const PropertyMap &map) const 
{
    //honestly, this is const, but the interface in FFMolecule
    //requires non-const...
    return InternalPotential::Molecule(molecule, 
                                       const_cast<InternalPotential&>(*this), 
                                       map);
}

/** Return the fully parameterised set of molecules from 'molecules', 
    using the supplied property map to find the properties containing the 
    required parameters.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InternalPotential::Molecules 
InternalPotential::parameterise(const MoleculeGroup &molecules,
                                const PropertyMap &map) const
{
    //honestly, this is const, but the interface in FFMolecules
    //requires non-const...
    return InternalPotential::Molecules(molecules, 
                                        const_cast<InternalPotential&>(*this), 
                                        map);
}

/** Return the coordinates of the atom 'atom' using the coordinates in 'cgroup_array' */
static const Vector& getCoords(const CGAtomIdx &atom,
                               const CoordGroup *cgroup_array)
{
    return cgroup_array[ atom.cutGroup() ].constData()[ atom.atom() ];
}

/** Calculate the energy caused by the physical terms (bond, angle, dihedral) */
void InternalPotential::calculatePhysicalEnergy(
                                         const GroupInternalParameters &group_params,
                                         const CoordGroup *cgroup_array,
                                         InternalPotential::Energy &energy,
                                         double scale_energy) const
{
    if (not group_params.bondPotential().isEmpty())
    {
        //bonds only use the symbol 'r', representing the bond length
        Values vals;
        const Symbol &r = InternalPotential::symbols().bond().r();
    
        int nbonds = group_params.bondPotential().count();
        const TwoAtomFunction *bonds_array = group_params.bondPotential().constData();
        
        double bndnrg = 0;
        
        for (int i=0; i<nbonds; ++i)
        {
            const TwoAtomFunction &bond = bonds_array[i];
            
            //get the coordinates of the two atoms in the bond
            const Vector &atom0 = getCoords(bond.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(bond.atom1(), cgroup_array);
            
            //calculate the interatomic distance, r
            double dist = Vector::distance(atom0, atom1);
            
            //substitute this distance into 'vals'
            vals.set( r, dist );
            
            //calculate the energy
            bndnrg += bond.function().evaluate(vals);
        }
        
        energy += BondEnergy( scale_energy * bndnrg );
    }
    
    if (not group_params.anglePotential().isEmpty())
    {
        //angles use the symbol 'theta', representing the angle size
        Values vals;
        const Symbol &theta = InternalPotential::symbols().angle().theta();
        
        int nangles = group_params.anglePotential().count();
        const ThreeAtomFunction *angles_array 
                                    = group_params.anglePotential().constData();
         
        double angnrg = 0;
                                                              
        for (int i=0; i<nangles; ++i)
        {
            const ThreeAtomFunction &angle = angles_array[i];
            
            Angle ang = Vector::angle( getCoords(angle.atom0(), cgroup_array),
                                       getCoords(angle.atom1(), cgroup_array),
                                       getCoords(angle.atom2(), cgroup_array) );
                           
            vals.set(theta, ang.to(radians));
            
            angnrg += angle.function().evaluate(vals);
        }
        
        energy += AngleEnergy( scale_energy * angnrg );
    }
    
    if (not group_params.dihedralPotential().isEmpty())
    {
        //angles use the symbol 'phi', representing the torsion angle
        Values vals;
        const Symbol &phi = InternalPotential::symbols().dihedral().phi();
        
        int ndihedrals = group_params.dihedralPotential().count();
        const FourAtomFunction *dihedrals_array 
                                    = group_params.dihedralPotential().constData();
         
        double dihnrg = 0;
                                                              
        for (int i=0; i<ndihedrals; ++i)
        {
            const FourAtomFunction &dihedral = dihedrals_array[i];
            
            Angle ang = Vector::dihedral( getCoords(dihedral.atom0(), cgroup_array),
                                          getCoords(dihedral.atom1(), cgroup_array),
                                          getCoords(dihedral.atom2(), cgroup_array),
                                          getCoords(dihedral.atom3(), cgroup_array) );
                           
            vals.set(phi, ang.to(radians));

            dihnrg += dihedral.function().evaluate(vals);
        }
        
        energy += DihedralEnergy( scale_energy * dihnrg );
    }
}

/** Calculate the energy caused by the non-physical terms (improper, Urey-Bradley) */
void InternalPotential::calculateNonPhysicalEnergy(
                                         const GroupInternalParameters &group_params,
                                         const CoordGroup *cgroup_array,
                                         InternalPotential::Energy &energy,
                                         double scale_energy) const
{
    if (not group_params.improperPotential().isEmpty())
    {
        Values vals;
        const Symbol &theta = InternalPotential::symbols().improper().theta();
        const Symbol &phi = InternalPotential::symbols().improper().phi();
        
        int nimpropers = group_params.improperPotential().count();
        const FourAtomFunction *impropers_array 
                                    = group_params.improperPotential().constData();
    
        double impnrg = 0;
        
        for (int i=0; i<nimpropers; ++i)
        {
            const FourAtomFunction &improper = impropers_array[i];
            
            Torsion torsion( getCoords(improper.atom0(), cgroup_array),
                             getCoords(improper.atom1(), cgroup_array),
                             getCoords(improper.atom2(), cgroup_array),
                             getCoords(improper.atom3(), cgroup_array) );
                             
            Angle theta_angle = torsion.improperAngle();
            Angle phi_angle = torsion.angle();
            
            vals.set(theta, theta_angle);
            vals.set(phi, phi_angle);
            
            impnrg += improper.function().evaluate(vals);
        }
        
        energy += ImproperEnergy( scale_energy * impnrg );
    }
    
    if (not group_params.ureyBradleyPotential().isEmpty())
    {
        Values vals;
        const Symbol &r = InternalPotential::symbols().ureyBradley().r();
        
        int nubs = group_params.ureyBradleyPotential().count();
        const TwoAtomFunction *ub_array 
                                = group_params.ureyBradleyPotential().constData();
                                
        double ubnrg = 0;
        
        for (int i=0; i<nubs; ++i)
        {
            const TwoAtomFunction &ub = ub_array[i];
            
            double dist = Vector::distance( getCoords(ub.atom0(), cgroup_array),
                                            getCoords(ub.atom1(), cgroup_array) );
                                            
            vals.set( r, dist );
            
            ubnrg += ub.function().evaluate(vals);
        }
        
        energy += UreyBradleyEnergy( scale_energy * ubnrg );
    }
}

/** Calculate the energy caused by the cross terms (improper, Urey-Bradley) */
void InternalPotential::calculateCrossEnergy(
                                         const GroupInternalParameters &group_params,
                                         const CoordGroup *cgroup_array,
                                         InternalPotential::Energy &energy,
                                         double scale_energy) const
{
    if (not group_params.stretchStretchPotential().isEmpty())
    {
        Values vals;
        const Symbol &r01 = InternalPotential::symbols().stretchStretch().r01();
        const Symbol &r21 = InternalPotential::symbols().stretchStretch().r21();
        
        int nss = group_params.stretchStretchPotential().count();
        const ThreeAtomFunction *ss_array 
                              = group_params.stretchStretchPotential().constData();
                              
        double ssnrg = 0;
        
        for (int i=0; i<nss; ++i)
        {
            const ThreeAtomFunction &ss = ss_array[i];
            
            const Vector &atom0 = getCoords(ss.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(ss.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(ss.atom2(), cgroup_array);
                           
            vals.set( r01, Vector::distance(atom0, atom1) );
            vals.set( r21, Vector::distance(atom2, atom1) );
            
            ssnrg += ss.function().evaluate(vals);
        }
        
        energy += StretchStretchEnergy( scale_energy * ssnrg );
    }

    if (not group_params.stretchBendPotential().isEmpty())
    {
        Values vals;
        const Symbol &r01 = InternalPotential::symbols().stretchBend().r01();
        const Symbol &r21 = InternalPotential::symbols().stretchBend().r21();
        const Symbol &theta = InternalPotential::symbols().stretchBend().theta();
        
        int nsb = group_params.stretchBendPotential().count();
        const ThreeAtomFunction *sb_array 
                              = group_params.stretchBendPotential().constData();
                              
        double sbnrg = 0;
        
        for (int i=0; i<nsb; ++i)
        {
            const ThreeAtomFunction &sb = sb_array[i];
            
            const Vector &atom0 = getCoords(sb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sb.atom2(), cgroup_array);
            
            vals.set( r01, Vector::distance(atom0, atom1) );
            vals.set( r21, Vector::distance(atom2, atom1) );
                
            vals.set( theta, Vector::angle(atom0, atom1, atom2) );
            
            sbnrg += sb.function().evaluate(vals);
        }
        
        energy += StretchBendEnergy( scale_energy * sbnrg );
    }

    if (not group_params.bendBendPotential().isEmpty())
    {
        Values vals;
        const Symbol &theta012 = InternalPotential::symbols().bendBend().theta012();
        const Symbol &theta213 = InternalPotential::symbols().bendBend().theta213();
        const Symbol &theta310 = InternalPotential::symbols().bendBend().theta310();
        
        int nbb = group_params.bendBendPotential().count();
        const FourAtomFunction *bb_array 
                              = group_params.bendBendPotential().constData();
                              
        double bbnrg = 0;
        
        for (int i=0; i<nbb; ++i)
        {
            const FourAtomFunction &bb = bb_array[i];
            
            const Vector &atom0 = getCoords(bb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(bb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(bb.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(bb.atom3(), cgroup_array);

            Vector v10 = atom0 - atom1;
            Vector v12 = atom2 - atom1;
            Vector v13 = atom3 - atom1;
        
            vals.set( theta012, Vector::angle(v12,v10) );
            vals.set( theta213, Vector::angle(v13,v12) );
            vals.set( theta310, Vector::angle(v10,v13) );
            
            bbnrg += bb.function().evaluate(vals);
        }
        
        energy += BendBendEnergy( scale_energy * bbnrg );
    }

    if (not group_params.stretchBendTorsionPotential().isEmpty())
    {
        Values vals;
    
        const Symbol &phi = InternalPotential::symbols().stretchBendTorsion().phi();
        const Symbol &r01 = InternalPotential::symbols().stretchBendTorsion().r01();
        const Symbol &r12 = InternalPotential::symbols().stretchBendTorsion().r12();
        const Symbol &r32 = InternalPotential::symbols().stretchBendTorsion().r32();
        const Symbol &r03 = InternalPotential::symbols().stretchBendTorsion().r03();
        const Symbol &theta012 
                = InternalPotential::symbols().stretchBendTorsion().theta012();
        const Symbol &theta321
                = InternalPotential::symbols().stretchBendTorsion().theta321();
        
        int nsbt = group_params.stretchBendTorsionPotential().count();
        const FourAtomFunction *sbt_array 
                          = group_params.stretchBendTorsionPotential().constData();
                              
        double sbtnrg = 0;
        
        for (int i=0; i<nsbt; ++i)
        {
            const FourAtomFunction &sbt = sbt_array[i];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);
            
            vals.set( phi, Vector::dihedral(atom0, atom1, atom2, atom3) );
            
            vals.set( r01, Vector::distance(atom0, atom1) );
            vals.set( r12, Vector::distance(atom1, atom2) );
            vals.set( r32, Vector::distance(atom2, atom3) );
            vals.set( r03, Vector::distance(atom0, atom3) );
            
            vals.set( theta012, Vector::angle(atom0, atom1, atom2) );
            vals.set( theta321, Vector::angle(atom3, atom2, atom1) );   
            
            sbtnrg += sbt.function().evaluate(vals);
        }
        
        energy += StretchBendTorsionEnergy( scale_energy * sbtnrg );
    }
}

/** Calculate the total intramolecular energy of 'molecule' and
    add it onto 'energy', optionally mutiplying it by a scale factor */
void InternalPotential::calculateEnergy(const InternalPotential::Molecule &molecule,
                                        InternalPotential::Energy &energy,
                                        double scale_energy) const
{
    if (scale_energy == 0)
        return;
        
    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];
        
        if (group_params.hasPhysicalParameters())
            this->calculatePhysicalEnergy(group_params, cgroup_array, 
                                          energy, scale_energy);
            
        if (group_params.hasNonPhysicalParameters())
            this->calculateNonPhysicalEnergy(group_params, cgroup_array,
                                             energy, scale_energy);
        
        if (group_params.hasCrossTerms())
            this->calculateCrossEnergy(group_params, cgroup_array,
                                       energy, scale_energy);
    }
}

static void addForce(const Vector &force, const CGAtomIdx &atom,
                     MolForceTable &forces)
{
    int idx = forces.map(atom.cutGroup());
    
    if (idx >= 0)
    {
        forces.data(idx)[atom.atom()] += force;
    }
}

/** Calculate the total bond force acting on the molecule 'molecule', and add it
    to the forces in 'forces', optionally scaled by 'scale_force' */
void InternalPotential::calculateBondForce(const InternalPotential::Molecule &molecule,
                                           MolForceTable &forces,
                                           double scale_force) const
{
    if (not molecule.parameters().hasBondParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &r = InternalPotential::symbols().bond().r();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any bonds?
        if ( group_params.bondForces().isEmpty() or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1())) )
        {
            continue;
        }
                      
        int nbonds = group_params.bondForces().count();
        const TwoAtomFunction *bonds_array = group_params.bondForces().constData();
        
        for (int j=0; j<nbonds; ++j)
        {
            const TwoAtomFunction &bond = bonds_array[j];
            
            const Vector &atom0 = getCoords(bond.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(bond.atom1(), cgroup_array);
            
            Vector v01 = atom0 - atom1;
            double dist = v01.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;
            
            v01 /= dist;
            vals.set( r, dist );

            //evaluate the force vector
            Vector force = -scale_force * bond.function().evaluate(vals) * v01;

            //add the force onto the forces array
            addForce(-force, bond.atom0(), forces);
            addForce(force, bond.atom1(), forces);
        }
    }
}

/** Calculate d theta / d r_i, d theta / d r_j and d theta / d r_k for 
    the passed r_i, r_j and r_k  (three points that make up an angle)
*/
static void d_theta_by_dr(const Vector &ri, const Vector &rj,
                          const Vector &rk,
                          
                          Vector &dtheta_by_dri, Vector &dtheta_by_drj,
                          Vector &dtheta_by_drk, 
                          
                          Angle &theta, double &dist_ij, double &dist_kj)
{
    Vector A = ri - rj;
    Vector B = rk - rj;
    
    double rA = A.length();
    double rB = B.length();
    
    dist_ij = rA;
    dist_kj = rB;
    
    if (rA * rB == 0)
    {
        //one of the bonds has zero length - cannot calculate angle force
        dtheta_by_dri = Vector(0);
        dtheta_by_drj = Vector(0);
        dtheta_by_drk = Vector(0);
        
        theta = Angle(0);
        return;
    }

    rA = 1 / rA;
    rB = 1 / rB;
    
    //work with normalised vectors to prevent numerical error
    A *= rA;
    B *= rB;
    
    const double A_B = Vector::dot(A,B);
    
    if (A_B == 0)
    {
        //bonds are parallel - cannot calculate an angle force
        dtheta_by_dri = Vector(0);
        dtheta_by_drj = Vector(0);
        dtheta_by_drk = Vector(0);
        
        theta = Angle(0);
        return;
    }

    dtheta_by_dri = rA * (B + A_B*A);
    dtheta_by_drk = rB * (A + A_B*B);
    
    dtheta_by_drj = -(dtheta_by_dri+dtheta_by_drk);
    
    theta = Angle( std::acos(A_B) );

    return;
}
                   
/** Calculate the total angle force acting on the molecule 'molecule', and add it
    to the forces in 'forces', optionally scaled by 'scale_force' 
    
*/
void InternalPotential::calculateAngleForce(const InternalPotential::Molecule &molecule,
                                            MolForceTable &forces,
                                            double scale_force) const
{
    if (not molecule.parameters().hasAngleParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &theta = InternalPotential::symbols().angle().theta();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any angles?
        if ( group_params.angleForces().isEmpty() or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2())) )
        {
            continue;
        }
                      
        int nangles = group_params.angleForces().count();
        const ThreeAtomFunction *angles_array = group_params.angleForces().constData();
        
        for (int j=0; j<nangles; ++j)
        {
            const ThreeAtomFunction &angle = angles_array[j];
            
            const Vector &atom0 = getCoords(angle.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(angle.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(angle.atom2(), cgroup_array);

            //-d V(theta) / dr  = -d V(theta) / dtheta  *  dtheta / dr
            
            //so first calculate d theta / dr
            Vector dtheta_by_d0, dtheta_by_d1, dtheta_by_d2;
            Angle t;
            double r01, r21;
            
            d_theta_by_dr(atom0, atom1, atom2,
                          dtheta_by_d0, dtheta_by_d1, dtheta_by_d2,
                          t, r01, r21);
                          
            //now calcualte -d V(theta) / d theta
            vals.set(theta, t);
            const double dv_by_dtheta = -scale_force * angle.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d0, angle.atom0(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, angle.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d2, angle.atom2(), forces);
        }
    }
}

static Vector cross(const Vector &v0, const Vector &v1)
{
    return Vector( v0.y()*v1.z() - v0.z()*v1.y(),
                   v0.z()*v1.x() - v0.x()*v1.z(),
                   v0.x()*v1.y() - v0.y()*v1.x() );
}

/** Calculate d phi / d r_i, d phi / d r_j, d phi / d r_k and d phi / d r_l.

    This uses the algorithm presented in Blondel and Karplus, J. Comp. Chem.
    Vol. 17, No. 9, 1132-1141, 1996, which shows how to calculate the derivative
    avoiding the singularity at phi = 0 (where sin(phi) = 0, so 1 / sin(phi) is
    undefined).
*/
static void d_phi_by_dr(const Vector &ri, const Vector &rj,
                        const Vector &rk, const Vector &rl,
                        
                        Vector &dphi_by_dri, Vector &dphi_by_drj,
                        Vector &dphi_by_drk, Vector &dphi_by_drl)
{
    Vector F = ri - rj;
    Vector G = rj - rk;
    Vector H = rl - rk;

    double rF = F.length();
    double rG = G.length();
    double rH = H.length();
    
    if (rF * rG * rH == 0)
    {
        //one of the bonds has zero length - cannot calculate
        //the dihedral force
        dphi_by_dri = Vector(0);
        dphi_by_drj = Vector(0);
        dphi_by_drk = Vector(0);
        dphi_by_drl = Vector(0);
        return;
    }
    
    //work with normalised vectors to prevent build-up of numerical error
    rF = 1 / rF;
    rG = 1 / rG;
    rH = 1 / rH;
    
    F *= rF;
    G *= rG;
    H *= rH;
    
    const double F_G = Vector::dot(F,G);
    const double H_G = Vector::dot(H,G);
    
    Vector A = ::cross(F, G);
    Vector B = ::cross(H, G);
    
    const double rA = A.length();
    const double rB = B.length();
    
    if (rA * rB == 0)
    {
        //atoms lie in a plane - cannot calculate the dihedral
        dphi_by_dri = Vector(0);
        dphi_by_drj = Vector(0);
        dphi_by_drk = Vector(0);
        dphi_by_drl = Vector(0);
        return;
    }

    A *= (1/rA);
    B *= (1/rB);
    
    const Vector dphi_by_dF = -rF * A;
    const Vector dphi_by_dH = rH * B;
    
    const Vector dphi_by_dG = rG * (F_G * A - H_G * B);

    dphi_by_dri = dphi_by_dF;
    dphi_by_drj = dphi_by_dG - dphi_by_dF;
    dphi_by_drk = -dphi_by_dG - dphi_by_dH;
    dphi_by_drl = dphi_by_dH;
}

/** Calculate d theta / d r_i, d theta / d r_j, d theta / d r_k and d theta / d r_l.

                                                  rk
    Four points, ri, rj, rk and rl               /
                                          ri--rj
                                                 \
                                                  rl
                                                  
    The theta angle is that between the vector ri-rj and the
    vector perpendicular to the plane rj, rk, rl
*/
static void d_theta_by_dr(const Vector &ri, const Vector &rj,
                          const Vector &rk, const Vector &rl,
                          Vector &dtheta_by_dri, Vector &dtheta_by_drj,
                          Vector &dtheta_by_drk, Vector &dtheta_by_drl)
{
    throw SireError::incomplete_code( QObject::tr(
        "Haven't yet implemented improper theta forces!"),
            CODELOC );
}

/** Calculate the total dihedral force acting on the molecule 'molecule', and add it
    to the forces in 'forces', optionally scaled by 'scale_force' 
    
*/
void InternalPotential::calculateDihedralForce(
                                    const InternalPotential::Molecule &molecule,
                                    MolForceTable &forces,
                                    double scale_force) const
{
    if (not molecule.parameters().hasDihedralParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &phi = InternalPotential::symbols().dihedral().phi();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any angles?
        if ( group_params.dihedralForces().isEmpty() or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2()) or
                  forces.selected(group_params.cgIdx3()) ) )
        {
            continue;
        }
                      
        int ndihedrals = group_params.dihedralForces().count();
        const FourAtomFunction *dihedrals_array 
                                    = group_params.dihedralForces().constData();
        
        for (int j=0; j<ndihedrals; ++j)
        {
            const FourAtomFunction &dihedral = dihedrals_array[j];
            
            const Vector &atom0 = getCoords(dihedral.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(dihedral.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(dihedral.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(dihedral.atom3(), cgroup_array);

            //-d V(phi) / dr  = -d V(phi) / dphi  *  dphi / dr
            
            //so first calculate d phi / dr
            Vector dphi_by_d0, dphi_by_d1, dphi_by_d2, dphi_by_d3;
            
            d_phi_by_dr(atom0, atom1, atom2, atom3,
                        dphi_by_d0, dphi_by_d1, dphi_by_d2, dphi_by_d3);
                          
            //now calcualte -d V(phi) / d phi
            vals.set(phi, Vector::dihedral(atom0,atom1,atom2,atom3).to(radians));
            double dv_by_dphi = -scale_force * dihedral.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dphi * dphi_by_d0, dihedral.atom0(), forces);
            addForce(dv_by_dphi * dphi_by_d1, dihedral.atom1(), forces);
            addForce(dv_by_dphi * dphi_by_d2, dihedral.atom2(), forces);
            addForce(dv_by_dphi * dphi_by_d3, dihedral.atom3(), forces);
        }
    }
}

/** Calculate the total improper force acting on the molecule 'molecule', and add it
    to the forces in 'forces', optionally scaled by 'scale_force' 
    
*/
void InternalPotential::calculateImproperForce(
                                    const InternalPotential::Molecule &molecule,
                                    MolForceTable &forces,
                                    double scale_force) const
{
    if (not molecule.parameters().hasImproperParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &phi = InternalPotential::symbols().improper().phi();
    const Symbol &theta = InternalPotential::symbols().improper().theta();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any angles?
        if ( (group_params.improper_Theta_Forces().isEmpty() and
              group_params.improper_Phi_Forces().isEmpty()) or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2()) or
                  forces.selected(group_params.cgIdx3()) ) )
        {
            continue;
        }
                      
        //do the theta forces first...
        int nimpropers = group_params.improper_Theta_Forces().count();
        const FourAtomFunction *impropers_array 
                                 = group_params.improper_Theta_Forces().constData();
        
        for (int j=0; j<nimpropers; ++j)
        {
            const FourAtomFunction &improper = impropers_array[j];
            
            const Vector &atom0 = getCoords(improper.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(improper.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(improper.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(improper.atom3(), cgroup_array);

            //d V(theta) / dr  = d V(theta) / dtheta  *  dtheta / dr
            
            //so first calculate d theta / dr
            Vector dtheta_by_d0, dtheta_by_d1, dtheta_by_d2, dtheta_by_d3;
            
            d_theta_by_dr(atom0, atom1, atom2, atom3,
                          dtheta_by_d0, dtheta_by_d1, dtheta_by_d2, dtheta_by_d3);
                          
            //now calculate d V(phi,theta) / d theta
            Torsion torsion(atom0, atom1, atom2, atom3);
                             
            vals.set(theta, torsion.improperAngle());
            vals.set(phi, torsion.angle());

            double dv_by_dtheta = scale_force * improper.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d0, improper.atom0(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, improper.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d2, improper.atom2(), forces);
            addForce(dv_by_dtheta * dtheta_by_d3, improper.atom3(), forces);
        }
        
        //now do the phi forces
        nimpropers = group_params.improper_Phi_Forces().count();
        impropers_array = group_params.improper_Phi_Forces().constData();
        
        for (int j=0; j<nimpropers; ++j)
        {
            const FourAtomFunction &improper = impropers_array[j];
            
            const Vector &atom0 = getCoords(improper.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(improper.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(improper.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(improper.atom3(), cgroup_array);

            //d V(phi) / dr  = d V(phi) / dphi  *  dphi / dr
            
            //so first calculate d phi / dr
            Vector dphi_by_d0, dphi_by_d1, dphi_by_d2, dphi_by_d3;
            
            d_phi_by_dr(atom0, atom1, atom2, atom3,
                        dphi_by_d0, dphi_by_d1, dphi_by_d2, dphi_by_d3);
                          
            //now calculate d V(phi,theta) / d phi
            Torsion torsion(atom0, atom1, atom2, atom3);
                             
            vals.set(theta, torsion.improperAngle());
            vals.set(phi, torsion.angle());

            double dv_by_dphi = scale_force * improper.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dphi * dphi_by_d0, improper.atom0(), forces);
            addForce(dv_by_dphi * dphi_by_d1, improper.atom1(), forces);
            addForce(dv_by_dphi * dphi_by_d2, improper.atom2(), forces);
            addForce(dv_by_dphi * dphi_by_d3, improper.atom3(), forces);
        }
    }
}
                            
/** Calculate the total Urey-Bradley force acting on the molecule 'molecule', 
    and add it to the forces in 'forces', optionally scaled by 'scale_force' */
void InternalPotential::calculateUBForce(const InternalPotential::Molecule &molecule,
                                         MolForceTable &forces,
                                         double scale_force) const
{
    if (not molecule.parameters().hasUreyBradleyParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &r = InternalPotential::symbols().ureyBradley().r();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any Urey Bradleys?
        if ( group_params.ureyBradleyForces().isEmpty() or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1())) )
        {
            continue;
        }
                      
        int nubs = group_params.ureyBradleyForces().count();
        const TwoAtomFunction *ubs_array = group_params.ureyBradleyForces().constData();
        
        for (int j=0; j<nubs; ++j)
        {
            const TwoAtomFunction &ub = ubs_array[j];
            
            const Vector &atom0 = getCoords(ub.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(ub.atom1(), cgroup_array);
            
            Vector v01 = atom0 - atom1;
            double dist = v01.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;
            
            v01 /= dist;
            vals.set( r, dist );

            //evaluate the force vector
            Vector force = scale_force * ub.function().evaluate(vals) * v01;

            //add the force onto the forces array
            addForce(force, ub.atom0(), forces);
            addForce(-force, ub.atom1(), forces);
        }
    }
}

/** Calculate the total stretch-stretch force acting on the molecule 'molecule', 
    and add it to the forces in 'forces', optionally scaled by 'scale_force' */
void InternalPotential::calculateSSForce(const InternalPotential::Molecule &molecule,
                                         MolForceTable &forces,
                                         double scale_force) const
{
    if (not molecule.parameters().hasStretchStretchParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &r01 = InternalPotential::symbols().stretchStretch().r01();
    const Symbol &r21 = InternalPotential::symbols().stretchStretch().r21();

    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group 
        //have any stretch-stretch forces?
        if ( (group_params.stretchStretch_R01_Forces().isEmpty() and
              group_params.stretchStretch_R21_Forces().isEmpty()) or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2())) )
        {
            continue;
        }
                      
        //calculate the forces along r01
        int nss = group_params.stretchStretch_R01_Forces().count();
        const ThreeAtomFunction *ss_array 
                            = group_params.stretchStretch_R01_Forces().constData();
        
        for (int j=0; j<nss; ++j)
        {
            const ThreeAtomFunction &ss = ss_array[j];
            
            const Vector &atom0 = getCoords(ss.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(ss.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(ss.atom2(), cgroup_array);
            
            Vector v = atom0 - atom1;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            vals.set( r01, dist );
            vals.set( r21, Vector::distance(atom1, atom2) );

            //evaluate the force vector
            Vector force = scale_force * ss.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, ss.atom0(), forces);
            addForce(-force, ss.atom1(), forces);
        }
                      
        //now calculate the forces along r21
        nss = group_params.stretchStretch_R21_Forces().count();
        ss_array = group_params.stretchStretch_R21_Forces().constData();
        
        for (int j=0; j<nss; ++j)
        {
            const ThreeAtomFunction &ss = ss_array[j];
            
            const Vector &atom0 = getCoords(ss.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(ss.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(ss.atom2(), cgroup_array);
            
            Vector v = atom2 - atom1;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            vals.set( r21, dist );
            vals.set( r01, Vector::distance(atom1, atom0) );

            //evaluate the force vector
            Vector force = scale_force * ss.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, ss.atom2(), forces);
            addForce(-force, ss.atom1(), forces);
        }
    }
}

/** Calculate all of the stretch-bend forces acting on the passed molecule,
    adding the forces onto the passed force table, optionally scaled
    by 'scale_force' */
void InternalPotential::calculateSBForce(const InternalPotential::Molecule &molecule,
                                         MolForceTable &forces,
                                         double scale_force) const
{
    if (not molecule.parameters().hasStretchBendParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &theta = InternalPotential::symbols().stretchBend().theta();
    const Symbol &r01 = InternalPotential::symbols().stretchBend().r01();
    const Symbol &r21 = InternalPotential::symbols().stretchBend().r21();
    
    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any forces?
        if ( (group_params.stretchBend_Theta_Forces().isEmpty() and
              group_params.stretchBend_R01_Forces().isEmpty() and
              group_params.stretchBend_R21_Forces().isEmpty()) or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2())) )
        {
            continue;
        }
                      
        //do the angle forces first
        int nsbs = group_params.stretchBend_Theta_Forces().count();
        const ThreeAtomFunction *sbs_array
                            = group_params.stretchBend_Theta_Forces().constData();
        
        for (int j=0; j<nsbs; ++j)
        {
            const ThreeAtomFunction &sb = sbs_array[j];
            
            const Vector &atom0 = getCoords(sb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sb.atom2(), cgroup_array);

            //d V(theta) / dr  = d V(theta) / dtheta  *  dtheta / dr
            
            //so first calculate d theta / dr
            Vector dtheta_by_d0, dtheta_by_d1, dtheta_by_d2;
            Angle t;
            double dist01, dist21;
            
            d_theta_by_dr(atom0, atom1, atom2,
                          dtheta_by_d0, dtheta_by_d1, dtheta_by_d2,
                          t, dist01, dist21);
                          
            //now calcualte d V(theta) / d theta
            vals.set(theta, t);
            vals.set(r01, dist01);
            vals.set(r21, dist21);
            double dv_by_dtheta = scale_force * sb.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d0, sb.atom0(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, sb.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d2, sb.atom2(), forces);
        }
        
        //now the r01 forces
        nsbs = group_params.stretchBend_R01_Forces().count();
        sbs_array = group_params.stretchBend_R01_Forces().constData();
        
        for (int j=0; j<nsbs; ++j)
        {
            const ThreeAtomFunction &sb = sbs_array[j];
            
            const Vector &atom0 = getCoords(sb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sb.atom2(), cgroup_array);

            Vector v = atom0 - atom1;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            vals.set( theta, Vector::angle(atom0, atom1, atom2) );
            vals.set( r01, dist );
            vals.set( r21, Vector::distance(atom1, atom2) );

            //evaluate the force vector
            Vector force = scale_force * sb.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, sb.atom0(), forces);
            addForce(-force, sb.atom1(), forces);
        }

        //now the r21 forces
        nsbs = group_params.stretchBend_R21_Forces().count();
        sbs_array = group_params.stretchBend_R21_Forces().constData();
        
        for (int j=0; j<nsbs; ++j)
        {
            const ThreeAtomFunction &sb = sbs_array[j];
            
            const Vector &atom0 = getCoords(sb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sb.atom2(), cgroup_array);

            Vector v = atom2 - atom1;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            vals.set( theta, Vector::angle(atom0, atom1, atom2) );
            vals.set( r21, dist );
            vals.set( r01, Vector::distance(atom1, atom0) );

            //evaluate the force vector
            Vector force = scale_force * sb.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, sb.atom2(), forces);
            addForce(-force, sb.atom1(), forces);
        }
    }
}

/** Calculate all of the bend-bend forces acting on the passed molecule,
    adding the forces onto the passed force table, optionally scaled
    by 'scale_force' */
void InternalPotential::calculateBBForce(const InternalPotential::Molecule &molecule,
                                         MolForceTable &forces,
                                         double scale_force) const
{
    if (not molecule.parameters().hasBendBendParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &theta012 = InternalPotential::symbols().bendBend().theta012();
    const Symbol &theta213 = InternalPotential::symbols().bendBend().theta213();
    const Symbol &theta310 = InternalPotential::symbols().bendBend().theta310();
    
    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any forces?
        if ( (group_params.bendBend_Theta012_Forces().isEmpty() and
              group_params.bendBend_Theta213_Forces().isEmpty() and
              group_params.bendBend_Theta310_Forces().isEmpty()) or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2())) )
        {
            continue;
        }
                      
        //do the theta012 forces first
        int nbbs = group_params.bendBend_Theta012_Forces().count();
        const FourAtomFunction *bbs_array
                            = group_params.bendBend_Theta012_Forces().constData();
        
        for (int j=0; j<nbbs; ++j)
        {
            const FourAtomFunction &bb = bbs_array[j];
            
            const Vector &atom0 = getCoords(bb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(bb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(bb.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(bb.atom3(), cgroup_array);

            //d V(theta) / dr  = d V(theta) / dtheta  *  dtheta / dr
            
            //so first calculate d theta / dr
            Vector dtheta_by_d0, dtheta_by_d1, dtheta_by_d2;
            Angle t;
            double dist01, dist21;
            
            d_theta_by_dr(atom0, atom1, atom2,
                          dtheta_by_d0, dtheta_by_d1, dtheta_by_d2,
                          t, dist01, dist21);
                          
            //now calcualte d V(theta) / d theta
            vals.set(theta012, t);
            vals.set(theta213, Vector::angle(atom2, atom1, atom3));
            vals.set(theta310, Vector::angle(atom3, atom1, atom0));

            double dv_by_dtheta = scale_force * bb.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d0, bb.atom0(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, bb.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d2, bb.atom2(), forces);
        }

        //now the theta213 forces first
        nbbs = group_params.bendBend_Theta213_Forces().count();
        bbs_array = group_params.bendBend_Theta213_Forces().constData();
        
        for (int j=0; j<nbbs; ++j)
        {
            const FourAtomFunction &bb = bbs_array[j];
            
            const Vector &atom0 = getCoords(bb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(bb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(bb.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(bb.atom3(), cgroup_array);

            //d V(theta) / dr  = d V(theta) / dtheta  *  dtheta / dr
            
            //so first calculate d theta / dr
            Vector dtheta_by_d2, dtheta_by_d1, dtheta_by_d3;
            Angle t;
            double dist21, dist31;
            
            d_theta_by_dr(atom2, atom1, atom3,
                          dtheta_by_d2, dtheta_by_d1, dtheta_by_d3,
                          t, dist21, dist31);
                          
            //now calcualte d V(theta) / d theta
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta213, t);
            vals.set(theta310, Vector::angle(atom3, atom1, atom0));

            double dv_by_dtheta = scale_force * bb.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d2, bb.atom2(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, bb.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d3, bb.atom3(), forces);
        }

        //finally the theta310 forces
        nbbs = group_params.bendBend_Theta310_Forces().count();
        bbs_array = group_params.bendBend_Theta310_Forces().constData();
        
        for (int j=0; j<nbbs; ++j)
        {
            const FourAtomFunction &bb = bbs_array[j];
            
            const Vector &atom0 = getCoords(bb.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(bb.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(bb.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(bb.atom3(), cgroup_array);

            //d V(theta) / dr  = d V(theta) / dtheta  *  dtheta / dr
            
            //so first calculate d theta / dr
            Vector dtheta_by_d3, dtheta_by_d1, dtheta_by_d0;
            Angle t;
            double dist31, dist01;
            
            d_theta_by_dr(atom3, atom1, atom0,
                          dtheta_by_d3, dtheta_by_d1, dtheta_by_d0,
                          t, dist31, dist01);
                          
            //now calcualte d V(theta) / d theta
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta213, Vector::angle(atom2, atom1, atom3));
            vals.set(theta310, t);

            double dv_by_dtheta = scale_force * bb.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d3, bb.atom3(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, bb.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d0, bb.atom0(), forces);
        }
    }
}

/** Calculate all of the stretch-bend-torsion forces acting on the passed molecule,
    adding the forces onto the passed force table, optionally scaled
    by 'scale_force' */
void InternalPotential::calculateSBTForce(const InternalPotential::Molecule &molecule,
                                          MolForceTable &forces,
                                          double scale_force) const
{
    if (not molecule.parameters().hasStretchBendTorsionParameters() or
            scale_force == 0)
        return;

    //get the array of all parameters for each group...
    const GroupInternalParameters *params_array 
                             = molecule.parameters().groupParameters().constData();
    
    int ngroups = molecule.parameters().groupParameters().count();

    //get the array of CoordGroups
    const CoordGroup *cgroup_array 
                             = molecule.parameters().atomicCoordinates().constData();

    Values vals;
    const Symbol &phi = InternalPotential::symbols().stretchBendTorsion().phi();

    const Symbol &r01 = InternalPotential::symbols().stretchBendTorsion().r01();
    const Symbol &r12 = InternalPotential::symbols().stretchBendTorsion().r12();
    const Symbol &r32 = InternalPotential::symbols().stretchBendTorsion().r32();
    const Symbol &r03 = InternalPotential::symbols().stretchBendTorsion().r03();
    
    const Symbol &theta012 
                = InternalPotential::symbols().stretchBendTorsion().theta012();
    const Symbol &theta321 
                = InternalPotential::symbols().stretchBendTorsion().theta321();
    
    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group_params = params_array[i];

        //Is this group in the force table, and does this group have any forces?
        if ( (group_params.stretchBendTorsion_Phi_Forces().isEmpty() and
              group_params.stretchBendTorsion_R01_Forces().isEmpty() and
              group_params.stretchBendTorsion_R12_Forces().isEmpty() and
              group_params.stretchBendTorsion_R32_Forces().isEmpty() and
              group_params.stretchBendTorsion_R03_Forces().isEmpty() and
              group_params.stretchBendTorsion_Theta012_Forces().isEmpty() and
              group_params.stretchBendTorsion_Theta321_Forces().isEmpty()) or
             not (forces.selected(group_params.cgIdx0()) or
                  forces.selected(group_params.cgIdx1()) or
                  forces.selected(group_params.cgIdx2())) )
        {
            continue;
        }
                      
        //do the angle forces first
        int nsbts = group_params.stretchBendTorsion_Phi_Forces().count();
        const FourAtomFunction *sbts_array
                        = group_params.stretchBendTorsion_Phi_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            //d V(phi) / dr  = d V(phi) / dphi  *  dphi / dr
            
            //so first calculate d phi / dr
            Vector dphi_by_d0, dphi_by_d1, dphi_by_d2, dphi_by_d3;
            
            d_phi_by_dr(atom0, atom1, atom2, atom3,
                        dphi_by_d0, dphi_by_d1, dphi_by_d2, dphi_by_d3);
                          
            //now calculate d V(phi,theta) / d phi
            Torsion torsion(atom0, atom1, atom2, atom3);
                             
            vals.set(phi, torsion.angle());
            vals.set(r01, Vector::distance(atom0, atom1));
            vals.set(r12, Vector::distance(atom1, atom2));
            vals.set(r32, Vector::distance(atom3, atom2));
            vals.set(r03, Vector::distance(atom0, atom3));
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta321, Vector::angle(atom3, atom2, atom1));

            double dv_by_dphi = scale_force * sbt.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dphi * dphi_by_d0, sbt.atom0(), forces);
            addForce(dv_by_dphi * dphi_by_d1, sbt.atom1(), forces);
            addForce(dv_by_dphi * dphi_by_d2, sbt.atom2(), forces);
            addForce(dv_by_dphi * dphi_by_d3, sbt.atom3(), forces);
        }
        
        //now the r01 forces
        nsbts = group_params.stretchBendTorsion_R01_Forces().count();
        sbts_array = group_params.stretchBendTorsion_R01_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            Vector v = atom0 - atom1;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            Torsion torsion(atom0, atom1, atom2, atom3);
            
            vals.set(phi, torsion.angle());
            vals.set(r01, dist);
            vals.set(r12, Vector::distance(atom1, atom2));
            vals.set(r32, Vector::distance(atom3, atom2));
            vals.set(r03, Vector::distance(atom0, atom3));
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta321, Vector::angle(atom3, atom2, atom1));

            //evaluate the force vector
            Vector force = scale_force * sbt.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, sbt.atom0(), forces);
            addForce(-force, sbt.atom1(), forces);
        }
        
        //now the r12 forces
        nsbts = group_params.stretchBendTorsion_R12_Forces().count();
        sbts_array = group_params.stretchBendTorsion_R12_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            Vector v = atom1 - atom2;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            Torsion torsion(atom0, atom1, atom2, atom3);
            
            vals.set(phi, torsion.angle());
            vals.set(r01, Vector::distance(atom0, atom1));
            vals.set(r12, dist);
            vals.set(r32, Vector::distance(atom3, atom2));
            vals.set(r03, Vector::distance(atom0, atom3));
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta321, Vector::angle(atom3, atom2, atom1));

            //evaluate the force vector
            Vector force = scale_force * sbt.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, sbt.atom1(), forces);
            addForce(-force, sbt.atom2(), forces);
        }
        
        //now the r32 forces
        nsbts = group_params.stretchBendTorsion_R32_Forces().count();
        sbts_array = group_params.stretchBendTorsion_R32_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            Vector v = atom3 - atom2;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            Torsion torsion(atom0, atom1, atom2, atom3);
            
            vals.set(phi, torsion.angle());
            vals.set(r01, Vector::distance(atom0, atom1));
            vals.set(r12, Vector::distance(atom1, atom2));
            vals.set(r32, dist);
            vals.set(r03, Vector::distance(atom0, atom3));
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta321, Vector::angle(atom3, atom2, atom1));

            //evaluate the force vector
            Vector force = scale_force * sbt.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, sbt.atom3(), forces);
            addForce(-force, sbt.atom2(), forces);
        }
        
        //now the r03 forces
        nsbts = group_params.stretchBendTorsion_R03_Forces().count();
        sbts_array = group_params.stretchBendTorsion_R03_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            Vector v = atom0 - atom3;
            double dist = v.length();
            
            if (dist == 0)
                //cannot calculate forces of overlapping atoms!
                continue;

            v /= dist;
            
            Torsion torsion(atom0, atom1, atom2, atom3);
            
            vals.set(phi, torsion.angle());
            vals.set(r01, Vector::distance(atom0, atom1));
            vals.set(r12, Vector::distance(atom1, atom2));
            vals.set(r32, Vector::distance(atom3, atom2));
            vals.set(r03, dist);
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta321, Vector::angle(atom3, atom2, atom1));

            //evaluate the force vector
            Vector force = scale_force * sbt.function().evaluate(vals) * v;

            //add the force onto the forces array
            addForce(force, sbt.atom0(), forces);
            addForce(-force, sbt.atom3(), forces);
        }
        
        //now the theta012 forces
        nsbts = group_params.stretchBendTorsion_Theta012_Forces().count();
        sbts_array = group_params.stretchBendTorsion_Theta012_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            Vector dtheta_by_d0, dtheta_by_d1, dtheta_by_d2;
            Angle t;
            double dist01, dist21;
            
            d_theta_by_dr(atom0, atom1, atom2,
                          dtheta_by_d0, dtheta_by_d1, dtheta_by_d2,
                          t, dist01, dist21);
            
            Torsion torsion(atom0, atom1, atom2, atom3);
            
            vals.set(phi, torsion.angle());
            vals.set(r01, Vector::distance(atom0, atom1));
            vals.set(r12, Vector::distance(atom1, atom2));
            vals.set(r32, Vector::distance(atom3, atom2));
            vals.set(r03, Vector::distance(atom0, atom3));
            vals.set(theta012, t);
            vals.set(theta321, Vector::angle(atom3, atom2, atom1));

            //evaluate the force vector
            double dv_by_dtheta = scale_force * sbt.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d0, sbt.atom0(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, sbt.atom1(), forces);
            addForce(dv_by_dtheta * dtheta_by_d2, sbt.atom2(), forces);
        }

        //now the theta321 forces
        nsbts = group_params.stretchBendTorsion_Theta321_Forces().count();
        sbts_array = group_params.stretchBendTorsion_Theta321_Forces().constData();
        
        for (int j=0; j<nsbts; ++j)
        {
            const FourAtomFunction &sbt = sbts_array[j];
            
            const Vector &atom0 = getCoords(sbt.atom0(), cgroup_array);
            const Vector &atom1 = getCoords(sbt.atom1(), cgroup_array);
            const Vector &atom2 = getCoords(sbt.atom2(), cgroup_array);
            const Vector &atom3 = getCoords(sbt.atom3(), cgroup_array);

            Vector dtheta_by_d3, dtheta_by_d2, dtheta_by_d1;
            Angle t;
            double dist32, dist12;
            
            d_theta_by_dr(atom3, atom2, atom1,
                          dtheta_by_d3, dtheta_by_d2, dtheta_by_d1,
                          t, dist32, dist12);
            
            Torsion torsion(atom0, atom1, atom2, atom3);
            
            vals.set(phi, torsion.angle());
            vals.set(r01, Vector::distance(atom0, atom1));
            vals.set(r12, Vector::distance(atom1, atom2));
            vals.set(r32, Vector::distance(atom3, atom2));
            vals.set(r03, Vector::distance(atom0, atom3));
            vals.set(theta012, Vector::angle(atom0, atom1, atom2));
            vals.set(theta321, t);

            //evaluate the force vector
            double dv_by_dtheta = scale_force * sbt.function().evaluate(vals);

            //add the force onto the forces array
            addForce(dv_by_dtheta * dtheta_by_d3, sbt.atom3(), forces);
            addForce(dv_by_dtheta * dtheta_by_d2, sbt.atom2(), forces);
            addForce(dv_by_dtheta * dtheta_by_d1, sbt.atom1(), forces);
        }
    }
}

/** Calculate the total force acting on the molecule 'molecule', and add it
    to the forces in 'forces', optionally scaled by 'scale_force' */
void InternalPotential::calculateForce(const InternalPotential::Molecule &molecule,
                                       MolForceTable &forces,
                                       double scale_force) const
{
    if (scale_force == 0)
        return;
    
    if (molecule.parameters().hasPhysicalParameters())
    {    
        calculateBondForce(molecule, forces, scale_force);
        calculateAngleForce(molecule, forces, scale_force);
        calculateDihedralForce(molecule, forces, scale_force);
    }
    
    if (molecule.parameters().hasNonPhysicalParameters())
    {
        calculateImproperForce(molecule, forces, scale_force);
        calculateUBForce(molecule, forces, scale_force);
    }
    
    if (molecule.parameters().hasCrossTerms())
    {
        calculateSSForce(molecule, forces, scale_force);
        calculateSBForce(molecule, forces, scale_force);
        calculateBBForce(molecule, forces, scale_force);
        calculateSBTForce(molecule, forces, scale_force);
    }
}

/** Calculate the forces in the molecule caused by the component represented
    by the symbol 'symbol' acting on the molecule 'molecule' and add this
    onto the forces in 'forces', optionally scaled by 'scale_force' 
    
    \throw SireFF::missing_component
*/
void InternalPotential::calculateForce(const InternalPotential::Molecule &molecule,
                                       MolForceTable &forces,
                                       const Symbol &symbol,
                                       const Components &components,
                                       double scale_force) const
{
    if (symbol == components.total())
        calculateForce(molecule, forces, scale_force);

    else if (symbol == components.bond())
        calculateBondForce(molecule, forces, scale_force);

    else if (symbol == components.angle())
        calculateAngleForce(molecule, forces, scale_force);

    else if (symbol == components.dihedral())
        calculateDihedralForce(molecule, forces, scale_force);

    else if (symbol == components.improper())
        calculateImproperForce(molecule, forces, scale_force);

    else if (symbol == components.ureyBradley())
        calculateUBForce(molecule, forces, scale_force);

    else if (symbol == components.stretchStretch())
        calculateSSForce(molecule, forces, scale_force);

    else if (symbol == components.stretchBend())
        calculateSBForce(molecule, forces, scale_force);

    else if (symbol == components.bendBend())
        calculateBBForce(molecule, forces, scale_force);

    else if (symbol == components.stretchBendTorsion())
        calculateSBTForce(molecule, forces, scale_force);

    else
        throw SireFF::missing_component( QObject::tr(
            "There is no component %1 in this potential. Available "
            "components are %2.")
                .arg(symbol.toString(), Sire::toString(components.symbols())),
                    CODELOC );
}

////////
//////// Implementation of InternalFF
////////

static const RegisterMetaType<InternalFF> r_internalff;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const InternalFF &internalff)
{
    writeHeader(ds, r_internalff, 2);
    
    SharedDataStream sds(ds);
    
    sds << internalff.mols << internalff.changed_mols
        << internalff.props
        << internalff.cljgroups
        << internalff.propmaps
        << static_cast<const G1FF&>(internalff);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      InternalFF &internalff)
{
    VersionID v = readHeader(ds, r_internalff);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> internalff.mols >> internalff.changed_mols
            >> internalff.props
            >> internalff.cljgroups
            >> internalff.propmaps
            >> static_cast<G1FF&>(internalff);
        
        internalff._pvt_updateName();
        internalff.isstrict = internalff.props.property("strict")
                                              .asABoolean();
        
        internalff.calc_14_nrgs = internalff.props.property("calculate14")
                                                  .asABoolean();
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> internalff.mols >> internalff.changed_mols
            >> internalff.props
            >> static_cast<G1FF&>(internalff);
        
        internalff.cljgroups.clear();
        internalff.propmaps.clear();
        
        internalff._pvt_updateName();
        internalff.isstrict = internalff.props.property("strict")
                                              .asABoolean();
        
        internalff.calc_14_nrgs = false;
    }
    else
        throw version_error(v, "1", r_internalff, CODELOC);
    
    return ds;
}

/** Constructor */
InternalFF::InternalFF()
           : ConcreteProperty<InternalFF,G1FF>(),
             FF3D(), InternalPotential(), calc_14_nrgs(false)
{
    props.setProperty("strict", wrap(true));
    props.setProperty("calculate14", wrap(calc_14_nrgs));
    props.setProperty("combiningRules", wrap("arithmetic"));
}

/** Construct a named internal forcefield */
InternalFF::InternalFF(const QString &name)
           : ConcreteProperty<InternalFF,G1FF>(),
             FF3D(), InternalPotential(), calc_14_nrgs(false)
{
    G1FF::setName(name);
    props.setProperty("strict", wrap(true));
    props.setProperty("calculate14", wrap(calc_14_nrgs));
    props.setProperty("combiningRules", wrap("arithmetic"));
}

/** Copy constructor */
InternalFF::InternalFF(const InternalFF &other)
           : ConcreteProperty<InternalFF,G1FF>(other),
             FF3D(other), InternalPotential(other),
             mols(other.mols), changed_mols(other.changed_mols),
             cljgroups(other.cljgroups),
             propmaps(other.propmaps),
             ffcomponents(other.ffcomponents),
             props(other.props),
             calc_14_nrgs(other.calc_14_nrgs)
{}

/** Destructor */
InternalFF::~InternalFF()
{}

/** Copy assignment operator */
InternalFF& InternalFF::operator=(const InternalFF &other)
{
    if (this != &other)
    {
        G1FF::operator=(other);
        FF3D::operator=(other);
        InternalPotential::operator=(other);
        
        mols = other.mols;
        changed_mols = other.changed_mols;
        ffcomponents = other.ffcomponents;
        props = other.props;
        
        cljgroups = other.cljgroups;
        propmaps = other.propmaps;
        calc_14_nrgs = other.calc_14_nrgs;
    }
    
    return *this;
}

/** Comparison operator */
bool InternalFF::operator==(const InternalFF &other) const
{
    return G1FF::operator==(other) and cljgroups == other.cljgroups and
           propmaps == other.propmaps and calc_14_nrgs == other.calc_14_nrgs;
}

/** Comparison operator */
bool InternalFF::operator!=(const InternalFF &other) const
{
    return G1FF::operator!=(other);
}

/** Function used to perform the work of changing the name of this 
    forcefield - this renames the component symbols and the molecule group */
void InternalFF::_pvt_updateName()
{
    ffcomponents = Components( this->name() );
    G1FF::_pvt_updateName();
}

/** Are we recording changes? */
bool InternalFF::recordingChanges() const
{
    return not (G1FF::isDirty() and changed_mols.isEmpty());
}

/** Record the change described in 'change' */
void InternalFF::recordChange(const InternalFF::ChangedMolecule &change)
{
    if (change.isEmpty())
        return;

    MolNum molnum = change.number();
    
    if (changed_mols.contains(molnum))
    {
        ChangedMolecule &old_change = changed_mols[molnum];
        
        if (old_change.oldMolecule() == change.newMolecule())
        {
            //we have reverted the change!
            changed_mols.remove(molnum);
            
            if (changed_mols.isEmpty())
                //there are no changes left, so this forcefield 
                //must now be clean
                G1FF::setClean();
            
            return;
        }
        else
        {
            //this is yet another change
            changed_mols[molnum].change( change.newMolecule() );
        }
    }
    else
    {
        changed_mols.insert(molnum, change);
    }

    G1FF::setDirty();
}

/** Set whether or not this strictly include terms that
    involve *only* selected atoms. Otherwise this includes
    terms that involve at least one selected atom */
bool InternalFF::setStrict(bool strict)
{
    if (isstrict == strict)
        return false;
        
    //change - this requires completely reparameterising all
    //partial molecules...
    try
    {
        isstrict = strict;
    
        const ChunkedVector<InternalFF::Molecule> &mols_array = mols.moleculesByIndex();
        const PropertyMap *property_array = mols.parameterNamesByIndex().constData();

        int nmols = mols.count();
        
        InternalFF::Molecules new_mols;
        
        for (int i=0; i<nmols; ++i)
        {
            new_mols.add( mols_array[i].molecule(), 
                          property_array[i], *this, false );
        }
        
        mols = new_mols;
        this->mustNowRecalculateFromScratch();
        
        for (QHash<MolNum,CLJ14Group>::iterator it = cljgroups.begin();
             it != cljgroups.end();
             ++it)
        {
            it.value().setStrict(isstrict);
        }
    }
    catch(...)
    {
        isstrict = not isstrict;
        throw;
    }
    
    props.setProperty( "strict", wrap(isstrict) );
    
    return true;
}

/** Turn on the calculate of 1-4 nonbonded terms */
void InternalFF::enable14Calculation()
{
    if (not calc_14_nrgs)
    {
        calc_14_nrgs = true;
        bool is_strict = this->isStrict();
        CLJFunction::COMBINING_RULES combining_rules = this->combiningRules();

        const ChunkedVector<InternalPotential::Molecule> &mols_array = mols.moleculesByIndex();
    
        for (int i=0; i<mols.count(); ++i)
        {
            cljgroups.insert( mols_array[i].number(),
                              CLJ14Group(mols_array[i].molecule(), combining_rules, is_strict,
                                         propmaps.value(mols_array[i].number(),
                                                        PropertyMap())) );
        }
        
        props.setProperty("calculate14", wrap(true));
    }
}

/** Disable calculation of the 1-4 nonbonded terms */
void InternalFF::disable14Calculation()
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
void InternalFF::setArithmeticCombiningRules(bool on)
{
    if (on)
        this->setCombiningRules( CLJFunction::ARITHMETIC );
    else
        this->setCombiningRules( CLJFunction::GEOMETRIC );
}

/** Turn on or off use of geometric combining rules when calculating
    the 1-4 nonbonded energy */
void InternalFF::setGeometricCombiningRules(bool on)
{
    if (on)
        this->setCombiningRules( CLJFunction::GEOMETRIC );
    else
        this->setCombiningRules( CLJFunction::ARITHMETIC );
}

/** Return the type of combining rules used when calculating the 1-4 
    nonbonded energy */
CLJFunction::COMBINING_RULES InternalFF::combiningRules() const
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
bool InternalFF::setCombiningRules(CLJFunction::COMBINING_RULES rules)
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
bool InternalFF::usingArithmeticCombiningRules() const
{
    return combiningRules() == CLJFunction::ARITHMETIC;
}

/** Return whether or not geometric combining rules are used for the 1-4 
    nonbonded energy calculation */
bool InternalFF::usingGeometricCombiningRules() const
{
    return combiningRules() == CLJFunction::GEOMETRIC;
}

/** Turn on or off the calculation of the 1-4 terms, returning whether
    or not this changes the forcefield */
bool InternalFF::setUse14Calculation(bool on)
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
bool InternalFF::uses14Calculation() const
{
    return calc_14_nrgs;
}

/** Set the property 'name' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
bool InternalFF::setProperty(const QString &name, const Property &value)
{
    if (name == QLatin1String("strict"))
    {
        return this->setStrict( value.asABoolean() );
    }
    else if (name == QLatin1String("combiningRules"))
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
            "InternalFF does not have a property called \"%1\" that "
            "can be changed. Available properties are [ strict ].")
                .arg(name), CODELOC );
}

/** Return the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& InternalFF::property(const QString &name) const
{
    return props.property(name);
}

/** Return whether this forcefield contains the property called 'name' */
bool InternalFF::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return the values of all of the properties of this forcefield */
const Properties& InternalFF::properties() const
{
    return props;
}

/** Calculate the energies in molecules in the passed energy table 
    caused by this potential, and add them onto the energies already
    in the energy table (optionally scaled by 'scale_energy') */
void InternalFF::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "InternalFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the energies of molecules in the passed energies table 
    caused by the component of this potential represented by 
    'symbol', and add them onto the energies already
    in the energy table (optionally scaled by 'scale_energy') */
void InternalFF::energy(EnergyTable &energytable, const Symbol &symbol,
			double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "InternalFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the forces acting on molecules in the passed force table 
    caused by this potential, and add them onto the forces already
    in the force table (optionally scaled by 'scale_force') */
void InternalFF::force(ForceTable &forcetable, double scale_force)
{
    if (scale_force == 0)
        return;
        
    int nforcemols = forcetable.count();
    MolForceTable *forcetable_array = forcetable.data();
    
    const ChunkedVector<InternalFF::Molecule> &mols_array = mols.moleculesByIndex();

    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (mols.contains(molnum))
        {
            InternalPotential::calculateForce(mols_array[mols.indexOf(molnum)],
                                              moltable, scale_force);
        }
    }
}

/** Calculate the forces acting on molecules in the passed force table 
    caused by the component of this potential represented by 
    'symbol', and add them onto the forces already
    in the force table (optionally scaled by 'scale_force') */
void InternalFF::force(ForceTable &forcetable, const Symbol &symbol,
                       double scale_force)
{
    if (scale_force == 0)
        return;
        
    int nforcemols = forcetable.count();
    MolForceTable *forcetable_array = forcetable.data();
    
    const ChunkedVector<InternalFF::Molecule> &mols_array = mols.moleculesByIndex();

    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (mols.contains(molnum))
        {
            InternalPotential::calculateForce(mols_array[mols.indexOf(molnum)],
                                              moltable, symbol,
                                              components(), scale_force);
        }
    }
}

/** Set it that the forcefield must now be recalculate from scratch */
void InternalFF::mustNowRecalculateFromScratch()
{
    //record that the forcefield is dirty
    G1FF::setDirty();
    
    //clear any delta information
    changed_mols.clear();
    
    for (QHash<MolNum,CLJ14Group>::iterator it = cljgroups.begin();
         it != cljgroups.end();
         ++it)
    {
        it.value().mustNowRecalculateFromScratch();
    }
}

/** Internal function used to get a handle on the forcefield components */
const FFComponent& InternalFF::_pvt_components() const
{
    return ffcomponents;
}

/** Recalculate the energy of the current state of this forcefield. This
    will recalculate the energy using the quickest possible route, e.g.
    if will only recalculate the energies of molecules that have changed
    since the last evaluation */
void InternalFF::recalculateEnergy()
{
    if (changed_mols.isEmpty())
    {
        //nothing appears to have changed, so lets recalculate
        //everything from scratch
        int nmols = mols.count();
        const ChunkedVector<InternalPotential::Molecule> &mols_array 
                                            = mols.moleculesByIndex();
        
        Energy total_nrg;
        
        for (int i=0; i<nmols; ++i)
        {
            InternalPotential::calculateEnergy( mols_array[i], total_nrg );
        }
        
        if (calc_14_nrgs)
        {
            double cnrg = 0;
            double ljnrg = 0;
        
            for (QHash<MolNum,CLJ14Group>::iterator it = cljgroups.begin();
                 it != cljgroups.end();
                 ++it)
            {
                boost::tuple<double,double> nrg = it.value().energy();
                cnrg += nrg.get<0>();
                ljnrg += nrg.get<1>();
            }
            
            total_nrg += Intra14Energy(cnrg, ljnrg);
        }
        
        this->components().setEnergy(*this, total_nrg);
    }
    else
    {
        //only some of the molecules have changed - calculate the
        //change in energy
        Energy old_nrg;
        Energy new_nrg;

        if (calc_14_nrgs)
        {
            //need to set clean so that we can get the old CLJ energy
            this->setClean();
            double old_cnrg = FF::energy(this->components().intra14Coulomb());
            double old_ljnrg = FF::energy(this->components().intra14LJ());
            this->setDirty();
            
            old_nrg += Intra14Energy(old_cnrg, old_ljnrg);
        }
        
        for (QHash<MolNum,ChangedMolecule>::const_iterator 
                                        it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (it->changedAll())
            {
                //must recalculate the total change in energy
                InternalPotential::calculateEnergy( it->oldMolecule(), old_nrg );
                InternalPotential::calculateEnergy( it->newMolecule(), new_nrg );
            }
            else
            {
                //can just calculte the change in energy of the changed parts
                InternalPotential::calculateEnergy( it->oldParts(), old_nrg );
                InternalPotential::calculateEnergy( it->newParts(), new_nrg );
            }
        }
        
        if (calc_14_nrgs)
        {
            double cnrg(0);
            double ljnrg(0);

            for (QHash<MolNum,CLJ14Group>::iterator it = cljgroups.begin();
                 it != cljgroups.end();
                 ++it)
            {
                boost::tuple<double,double> nrg = it.value().energy();
                cnrg += nrg.get<0>();
                ljnrg += nrg.get<1>();
            }

            new_nrg += Intra14Energy(cnrg, ljnrg);
        }

        this->components().changeEnergy(*this, new_nrg - old_nrg);
        
        //can now forget about the changes :-)
        changed_mols.clear();
    }
    
    this->setClean();
}

/** Record that the molecule in 'mol' has been added, using the 
    passed property map to get the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void InternalFF::_pvt_added(const PartialMolecule &molecule, const PropertyMap &map)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols.add(molecule, map, *this, true);
        this->recordChange(mol);
    }
    else
    {
        mols.add(molecule, map, *this, false);
    }
    
    if (not map.isDefault())
    {
        propmaps.insert(molecule.number(), map);
    }
    
    if (calc_14_nrgs)
    {
        if (cljgroups.contains(molecule.number()))
        {
            cljgroups[molecule.number()].add(molecule);
        }
        else
        {
            cljgroups.insert(molecule.number(),
                             CLJ14Group(molecule, combiningRules(), isStrict(), map));
        }
    }
}

/** Record the fact that the molecule 'mol' has been removed from this forcefield */
void InternalFF::_pvt_removed(const PartialMolecule &molecule)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols.remove(molecule, *this, true);
        this->recordChange(mol);
    }
    else
    {
        mols.remove(molecule, *this, false);
    }
    
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
void InternalFF::_pvt_changed(const SireMol::Molecule &molecule, bool auto_update)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols.change(molecule, *this, true);
        this->recordChange(mol);
    }
    else
    {
        mols.change(molecule, *this, false);
    }
    
    if (calc_14_nrgs)
    {
        if (cljgroups.contains(molecule.number()))
        {
            cljgroups[molecule.number()].update(molecule);
        }
    }
}

/** Record that the provided list of molecules have changed 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void InternalFF::_pvt_changed(const QList<SireMol::Molecule> &molecules, bool auto_update)
{
    InternalFF::Molecules old_mols = mols;
    QHash<MolNum,ChangedMolecule> old_changed_mols = changed_mols;
    QHash<MolNum,CLJ14Group> old_cljgroups = cljgroups;

    try
    {
        if (this->recordingChanges())
        {   
            foreach (const SireMol::Molecule &molecule, molecules)
            {
                ChangedMolecule change = mols.change(molecule, *this, true);
                this->recordChange(change);
            }
        }
        else
        {
            foreach (const SireMol::Molecule &molecule, molecules)
            {
                mols.change(molecule, *this, false);
            }
        }
        
        if (calc_14_nrgs)
        {
            foreach (const SireMol::Molecule &molecule, molecules)
            {
                if (cljgroups.contains(molecule.number()))
                {
                    cljgroups[molecule.number()].update(molecule);
                }
            }
        }
    }
    catch(...)
    {
        mols = old_mols;
        changed_mols = old_changed_mols;
        cljgroups = old_cljgroups;
        throw;
    }
}

/** Record that all of the molecules have been removed */
void InternalFF::_pvt_removedAll()
{
    mols.clear();
    cljgroups.clear();
    this->mustNowRecalculateFromScratch();
}
    
/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */       
bool InternalFF::_pvt_wouldChangeProperties(MolNum molnum, 
                                            const PropertyMap &map) const
{
    if (mols.wouldChangeProperties(molnum, map))
        return true;

    else if (cljgroups.contains(molnum))
        return cljgroups[molnum].wouldChangeProperties(map);

    else
        return false;
}

const char* InternalFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalFF>() );
}

InternalFF* InternalFF::clone() const
{
    return new InternalFF(*this);
}

void InternalFF::field(FieldTable &fieldtable, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}

void InternalFF::field(FieldTable &fieldtable, const Symbol &component,
                       double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}
           
void InternalFF::potential(PotentialTable &potentialtable, double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}

void InternalFF::potential(PotentialTable &potentialtable, const Symbol &component,
                          double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}

void InternalFF::field(FieldTable &fieldtable, const Probe &probe, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}

void InternalFF::field(FieldTable &fieldtable, const Symbol &component,
                       const Probe &probe, double scale_field)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}

void InternalFF::potential(PotentialTable &potentialtable, const Probe &probe,
                           double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}

void InternalFF::potential(PotentialTable &potentialtable, const Symbol &component,
                           const Probe &probe, double scale_potential)
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the field of an InternalFF has yet to "
                "be implemented."), CODELOC );
}
