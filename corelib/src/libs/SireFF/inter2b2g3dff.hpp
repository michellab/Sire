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

#ifndef SIREFF_INTER2B2G3DFF_HPP
#define SIREFF_INTER2B2G3DFF_HPP

#include "ff3d.h"
#include "inter2b2gff.hpp"

SIRE_BEGIN_HEADER

namespace SireFF
{

/** This class provides an intermolecular non-bonded potential
    that can work with an two-body, three-dimensional potential
    (provided by the template class 'Potential') that
    calculates the intermolecular potential between
    two groups of molecules. This is a 3D specific forcefield, 
    so provides functions that can calculate the 3D forces acting on the atoms.
    
    @author Christopher Woods
*/
template<class Potential>
class Inter2B2G3DFF 
                : public SireBase::ConcreteProperty< Inter2B2G3DFF<Potential>,
                                                     Inter2B2GFF<Potential> >, 
                  public FF3D
{
public:
    Inter2B2G3DFF();
    Inter2B2G3DFF(const QString &name);
    
    Inter2B2G3DFF(const Inter2B2G3DFF<Potential> &other);
    
    ~Inter2B2G3DFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Inter2B2G3DFF<Potential>& operator=(const Inter2B2G3DFF<Potential> &other);
    
    bool operator==(const Inter2B2G3DFF<Potential> &other) const;
    bool operator!=(const Inter2B2G3DFF<Potential> &other) const;
    
    Inter2B2G3DFF<Potential>* clone() const;

    SireUnits::Dimension::MolarEnergy energy();
    SireUnits::Dimension::MolarEnergy energy(const Symbol &component);
    
    void energy(EnergyTable &energytable, double scale_energy=1);
    
    void energy(EnergyTable &energytable, const Symbol &symbol,
		double scale_energy=1);

    void force(ForceTable &forcetable, double scale_force=1);
    
    void force(ForceTable &forcetable, const Symbol &symbol,
               double scale_force=1);
               
    void field(FieldTable &fieldtable, double scale_field=1);
    
    void field(FieldTable &fieldtable, const Symbol &component,
               double scale_field=1);
               
    void potential(PotentialTable &potentialtable, double scale_potential=1);
    
    void potential(PotentialTable &potentialtable, const Symbol &component,
                   double scale_potential=1);

    void field(FieldTable &fieldtable, const Probe &probe, double scale_field=1);
    
    void field(FieldTable &fieldtable, const Symbol &component,
               const Probe &probe, double scale_field=1);
               
    void potential(PotentialTable &potentialtable, const Probe &probe,
                   double scale_potential=1);
    
    void potential(PotentialTable &potentialtable, const Symbol &component,
                   const Probe &probe, double scale_potential=1);

protected:
    typedef typename Inter2B2GFF<Potential>::Energy Energy;
    typedef typename Inter2B2GFF<Potential>::EnergyWorkspace EnergyWorkspace;

    typedef typename Inter2B2GFF<Potential>::Molecules Molecules;
    typedef typename Inter2B2GFF<Potential>::Molecule Molecule;

    typedef typename Inter2B2GFF<Potential>::ChangedMolecule ChangedMolecule;

    void recalculateEnergy();
};

#if !defined(SIRE_SKIP_INLINE_FUNCTIONS) \
    && !defined(InterCLJFF_hpp__pyplusplus_wrapper) \
    && !defined(InterCLJFFBase_hpp__pyplusplus_wrapper) \
    && !defined(CLJPotentialInterface_InterCLJPotential__hpp__pyplusplus_wrapper) \
    && !defined(InterGroupCLJFFBase_hpp__pyplusplus_wrapper) \
    && !defined(InterGroupCLJFF_hpp__pyplusplus_wrapper)
/** Constructor (without giving the forcefield a name!) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B2G3DFF<Potential>::Inter2B2G3DFF()
  : SireBase::ConcreteProperty< Inter2B2G3DFF<Potential>,Inter2B2GFF<Potential> >(), 
    FF3D()
{}

/** Construct, giving this forcefield a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B2G3DFF<Potential>::Inter2B2G3DFF(const QString &name)
  : SireBase::ConcreteProperty< Inter2B2G3DFF<Potential>,Inter2B2GFF<Potential> >(name), 
    FF3D()
{}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B2G3DFF<Potential>::Inter2B2G3DFF(const Inter2B2G3DFF<Potential> &other)
  : SireBase::ConcreteProperty< Inter2B2G3DFF<Potential>,Inter2B2GFF<Potential> >(other), 
    FF3D(other)
{}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B2G3DFF<Potential>::~Inter2B2G3DFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B2G3DFF<Potential>& 
Inter2B2G3DFF<Potential>::operator=(const Inter2B2G3DFF<Potential> &other)
{
    Inter2B2GFF<Potential>::operator=(other);
    FF3D::operator=(other);
    
    return *this;
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2B2G3DFF<Potential>::operator==(const Inter2B2G3DFF<Potential> &other) const
{
    return Inter2B2GFF<Potential>::operator==(other);
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2B2G3DFF<Potential>::operator!=(const Inter2B2G3DFF<Potential> &other) const
{
    return Inter2B2GFF<Potential>::operator!=(other);
}
    
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Inter2B2G3DFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Inter2B2G3DFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Inter2B2G3DFF<Potential>::what() const
{
    return Inter2B2G3DFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B2G3DFF<Potential>* Inter2B2G3DFF<Potential>::clone() const
{
    return new Inter2B2G3DFF<Potential>(*this);
}

/** Recalculate the energy of the current state of this forcefield. This
    will recalculate the energy using the quickest possible route, e.g.
    if will only recalculate the energies of molecules that have changed
    since the last evaluation */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::recalculateEnergy()
{
    int nmols0 = this->mols[0].count();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();

    int nmols1 = this->mols[1].count();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();

    //tell the potential that we are starting an evaluation
    Potential::startEvaluation();

    try
    {

    const ChunkedVector<SireVol::AABox> &aaboxes0_array 
                                         = this->mols[0].aaBoxesByIndex();

    const ChunkedVector<SireVol::AABox> &aaboxes1_array 
                                         = this->mols[1].aaBoxesByIndex();

    const SireVol::Space &spce = this->space();
    const double cutoff = this->switchingFunction().cutoffDistance();

    if (this->changed_mols[0].isEmpty() and this->changed_mols[1].isEmpty())
    {
        //we are not recording changes, so we have to assume that
        //everything has changed. Recalculate the total energy from scratch
        EnergyWorkspace workspace;
        Energy total_nrg;

        //loop over all pairs of molecules
        for (int i=0; i<nmols0; ++i)
        {
            const typename Potential::Molecule &mol0 = mols0_array.at(i);
            const SireVol::AABox &aabox0 = aaboxes0_array.at(i);
        
            for (int j=0; j<nmols1; ++j)
            {
                if (not spce.beyond(cutoff, aabox0, aaboxes1_array.at(j)))
                {                
                    const typename Potential::Molecule &mol1 = mols1_array.at(j);
                    Potential::calculateEnergy(mol0, mol1, total_nrg, workspace);
                }
            }
        }
        
        //set the energy
        this->components().setEnergy(*this, total_nrg);
    }
    else
    {
        //just calculate the changes in energy
        Energy old_nrg;
        Energy new_nrg;

        EnergyWorkspace old_workspace;
        EnergyWorkspace new_workspace;
        
        if (this->changed_mols[0].count() == 1 and this->changed_mols[1].isEmpty())
        {
            //only need to calculate the energy of changed_mol0 with changed_mols1
            const typename Potential::ChangedMolecule &changed_mol = 
                                            *(this->changed_mols[0].constBegin());
                         
            const SireVol::AABox &new_box = changed_mol.newMolecule().aaBox();
            const SireVol::AABox &old_box = changed_mol.oldMolecule().aaBox();

            for (int i=0; i<nmols1; ++i)
            {
                const typename Potential::Molecule &mol = mols1_array.at(i);
                const SireVol::AABox &aabox = aaboxes1_array.at(i);
                
                if (not spce.beyond(cutoff, old_box, aabox))
                {
                    Potential::calculateEnergy(mol, changed_mol.oldParts(),
                                               old_nrg, old_workspace);
                }

                if (not spce.beyond(cutoff, new_box, aabox))
                {
                    Potential::calculateEnergy(mol, changed_mol.newParts(),
                                               new_nrg, new_workspace);
                }
            }
        }
        else if (this->changed_mols[0].isEmpty() and this->changed_mols[1].count() == 1)
        {
            //only need to calculate the energy of changed_mol0 with changed_mols1
            const typename Potential::ChangedMolecule &changed_mol = 
                                            *(this->changed_mols[1].constBegin());
                         
            const SireVol::AABox &new_box = changed_mol.newMolecule().aaBox();
            const SireVol::AABox &old_box = changed_mol.oldMolecule().aaBox();

            for (int i=0; i<nmols0; ++i)
            {
                const typename Potential::Molecule &mol = mols0_array.at(i);
                const SireVol::AABox &aabox = aaboxes0_array.at(i);
                
                if (not spce.beyond(cutoff, old_box, aabox))
                {
                    Potential::calculateEnergy(mol, changed_mol.oldParts(),
                                               old_nrg, old_workspace);
                }

                if (not spce.beyond(cutoff, new_box, aabox))
                {
                    Potential::calculateEnergy(mol, changed_mol.newParts(),
                                               new_nrg, new_workspace);
                }
            }
        }
        else 
        {
            if (not this->changed_mols[1].isEmpty())
            {
                for (int i=0; i<nmols0; ++i)
                {
                    const typename Potential::Molecule &mol0 = mols0_array.at(i);
                    
                    if (this->changed_mols[0].contains(mol0.number()))
                        //this molecule has changed as well - deal with 
                        //this later...
                        continue;
                        
                    for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                                  it = this->changed_mols[1].constBegin();
                         it != this->changed_mols[1].constEnd();
                         ++it)
                    {
                        Potential::calculateEnergy(mol0, it->oldParts(),
                                                   old_nrg, old_workspace);
                                                   
                        Potential::calculateEnergy(mol0, it->newParts(),
                                                   new_nrg, new_workspace);
                    }
                }
            }
            
            if (not this->changed_mols[0].isEmpty())
            {
                for (int i=0; i<nmols1; ++i)
                {
                    const typename Potential::Molecule &mol1 = mols1_array.at(i);
                    
                    if (this->changed_mols[1].contains(mol1.number()))
                        //this molecule has changed as well - deal with
                        //this later...
                        continue;
                        
                    for (typename QHash<MolNum,ChangedMolecule>::const_iterator it
                                                  = this->changed_mols[0].constBegin();
                         it != this->changed_mols[0].constEnd();
                         ++it)
                    {
                        Potential::calculateEnergy(mol1, it->oldParts(),
                                                   old_nrg, old_workspace);
                                                   
                        Potential::calculateEnergy(mol1, it->newParts(),
                                                   new_nrg, new_workspace);
                    }
                }
            }
            
            //now calculate the energy change between the changed molecules
            //of group0 and the changed molecules of group1
            if ( not (this->changed_mols[0].isEmpty() and 
                      this->changed_mols[1].isEmpty()) )
            {
                for (typename QHash<MolNum,ChangedMolecule>::const_iterator it0
                                                   = this->changed_mols[0].constBegin();
                     it0 != this->changed_mols[0].constEnd();
                     ++it0)
                {
                    bool changed_all_0 = it0->changedAll();
                
                    for (typename QHash<MolNum,ChangedMolecule>::const_iterator it1
                                                   = this->changed_mols[1].constBegin();
                         it1 != this->changed_mols[1].constEnd();
                         ++it1)
                    {
                        if (changed_all_0 or it1->changedAll())
                        {
                            //need to do the whole molecule calculation
                            Potential::calculateEnergy( it0->oldMolecule(),
                                                        it1->oldMolecule(),
                                                        old_nrg, old_workspace );
                                                        
                            Potential::calculateEnergy( it0->newMolecule(),
                                                        it1->newMolecule(),
                                                        new_nrg, new_workspace );
                        }
                        else
                        {
                            //part of mol0 has changed and part of mol1 has changed
                            // - calculate the energies of the changes with 
                            //   the whole of the other molecule...
                            Potential::calculateEnergy( it0->oldParts(),
                                                        it1->oldMolecule(),
                                                        old_nrg, old_workspace );
                                                        
                            Potential::calculateEnergy( it0->newParts(),
                                                        it1->newMolecule(),
                                                        new_nrg, new_workspace );
                                                        
                            Potential::calculateEnergy( it1->oldParts(),
                                                        it0->oldMolecule(),
                                                        old_nrg, old_workspace );
                                                        
                            Potential::calculateEnergy( it1->newParts(),
                                                        it0->newMolecule(),
                                                        new_nrg, new_workspace );
                                                        
                            //now we need to subtract the double counted changed
                            //parts of mol0 with changed parts of mol1
                            Potential::calculateEnergy( it0->oldParts(),
                                                        it1->oldParts(),
                                                        old_nrg, old_workspace, -1 );
                                                        
                            Potential::calculateEnergy( it0->newParts(),
                                                        it1->newParts(),
                                                        new_nrg, new_workspace, -1 );
                        }
                    }
                }
            }
        }
         
        //change the energy
        this->components().changeEnergy(*this, new_nrg - old_nrg);
        
        //clear the changed molecules
        this->changed_mols[0].clear();
        this->changed_mols[1].clear();
    }

    Potential::finishedEvaluation();
    
    this->setClean();
    
    }
    catch(...)
    {
        Potential::finishedEvaluation();
        throw;
    }
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
SireUnits::Dimension::MolarEnergy Inter2B2G3DFF<Potential>::energy()
{
    return Inter2B2GFF<Potential>::energy();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
SireUnits::Dimension::MolarEnergy Inter2B2G3DFF<Potential>::energy(const Symbol &component)
{
    return Inter2B2GFF<Potential>::energy(component);
}

/** Calculate the energies of the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the energies present
    in the energy table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "Inter2B2G3DFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the energies acting on the molecules in the passed energytable  
    caused by the component of this forcefield represented by 'symbol',
    adding this energy onto the existing energies in the forcetable (optionally
    multiplied by 'scale_energy' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::energy(EnergyTable &energytable, const Symbol &symbol, double scale_energy)
{
    if (scale_energy == 0)
        return;

    int nenergymols = energytable.count();
    int nmols0 = this->mols[0].count();
    int nmols1 = this->mols[1].count();
    
    typename Potential::EnergyWorkspace workspace;
    
    MolEnergyTable *energytable_array = energytable.data();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<nenergymols; ++i)
    {
        MolEnergyTable &moltable = energytable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum))
        {
            //calculate the energies on this molecule caused by group0
            int imol = this->mols[0].indexOf(molnum);
            const typename Potential::Molecule &mol0 = mols0_array[imol];
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                
		Potential::calculateEnergy(mol0, mol1, moltable, workspace, scale_energy);		
            }
        }
        
        if (this->mols[1].contains(molnum))
        {
            //calculate the energies on this molecule caused by group1
            int imol = this->mols[1].indexOf(molnum);
            const typename Potential::Molecule &mol0 = mols1_array[imol];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol1 = mols0_array[j];
                
 		Potential::calculateEnergy(mol0, mol1, moltable, workspace, scale_energy);
            }
        }
    }

}

/** Calculate the forces acting on the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the forces present
    in the force table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::force(ForceTable &forcetable, double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    int nmols0 = this->mols[0].count();
    int nmols1 = this->mols[1].count();
    
    typename Potential::ForceWorkspace workspace;
    
    MolForceTable *forcetable_array = forcetable.data();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum))
        {
            //calculate the forces on this molecule caused by group0
            int imol = this->mols[0].indexOf(molnum);
            const typename Potential::Molecule &mol0 = mols0_array[imol];
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                
                Potential::calculateForce(mol0, mol1, moltable,
                                          workspace, scale_force);
            }
        }
        
        if (this->mols[1].contains(molnum))
        {
            //calculate the forces on this molecule caused by group1
            int imol = this->mols[1].indexOf(molnum);
            const typename Potential::Molecule &mol0 = mols1_array[imol];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol1 = mols0_array[j];
                
                Potential::calculateForce(mol0, mol1, moltable,
                                          workspace, scale_force);
            }
        }
    }
}

/** Calculate the force acting on the molecules in the passed forcetable  
    caused by the component of this forcefield represented by 'symbol',
    adding this force onto the existing forces in the forcetable (optionally
    multiplied by 'scale_force' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::force(ForceTable &forcetable, const Symbol &symbol,
                                   double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    int nmols0 = this->mols[0].count();
    int nmols1 = this->mols[1].count();
    
    typename Potential::ForceWorkspace workspace;
    
    MolForceTable *forcetable_array = forcetable.data();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum))
        {
            //calculate the forces on this molecule caused by group0
            int imol = this->mols[0].indexOf(molnum);
            const typename Potential::Molecule &mol0 = mols0_array[imol];
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                
                Potential::calculateForce(mol0, mol1, moltable,
                                          symbol, this->components(),
                                          workspace, scale_force);
            }
        }
        
        if (this->mols[1].contains(molnum))
        {
            //calculate the forces on this molecule caused by group1
            int imol = this->mols[1].indexOf(molnum);
            const typename Potential::Molecule &mol0 = mols1_array[imol];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol1 = mols0_array[j];
                
                Potential::calculateForce(mol0, mol1, moltable,
                                          symbol, this->components(),
                                          workspace, scale_force);
            }
        }
    }
}

/** Calculate the fields acting at the points in the passed fieldtable
    that arise from this forcefield, and add them onto the fields present
    in the field table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::field(FieldTable &fieldtable, const Probe &prob,
                                     double scale_field)
{
    if (scale_field == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int ngrids = fieldtable.nGrids();
    const int nfieldmols = fieldtable.nMolecules();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    
    typename Potential::FieldWorkspace workspace;
    
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();

    if (nfieldmols > 0)
    {
        MolFieldTable *moltable_array = fieldtable.moleculeData();
    
        for (int i=0; i<nfieldmols; ++i)
        {
            MolFieldTable &moltable = moltable_array[i];
            
            MolNum molnum = moltable.molNum();
            
            if (this->mols[0].contains(molnum))
            {
                //calculate the field on this molecule caused by group0
                int imol = this->mols[0].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols0_array[imol];
                
                for (int j=0; j<nmols1; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols1_array[j];
                    
                    Potential::calculateField(mol0, mol1, probe, moltable,
                                              workspace, scale_field);
                }
            }
            
            if (this->mols[1].contains(molnum))
            {
                //calculate the fields on this molecule caused by group1
                int imol = this->mols[1].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols1_array[imol];
                
                for (int j=0; j<nmols0; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols0_array[j];
                    
                    Potential::calculateField(mol0, mol1, probe, moltable,
                                              workspace, scale_field);
                }
            }
        }
    }
    
    if (ngrids > 0)
    {
        GridFieldTable *gridtable_array = fieldtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridFieldTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol0 = mols0_array[j];
                    
                Potential::calculateField(mol0, probe, gridtable,
                                          workspace, scale_field);
            }

            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                    
                Potential::calculateField(mol1, probe, gridtable,
                                          workspace, scale_field);
            }
        }
    }
}

/** Calculate the field acting on the points in the passed forcetable  
    caused by the component of this forcefield represented by 'component',
    adding this field onto the existing fields in the forcetable (optionally
    multiplied by 'scale_field' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                     const Probe &prob, double scale_field)
{
    if (scale_field == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int nfieldmols = fieldtable.nMolecules();
    const int ngrids = fieldtable.nGrids();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    
    typename Potential::FieldWorkspace workspace;
    
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    if (nfieldmols > 0)
    {
        MolFieldTable *fieldtable_array = fieldtable.moleculeData();

        for (int i=0; i<nfieldmols; ++i)
        {
            MolFieldTable &moltable = fieldtable_array[i];
        
            MolNum molnum = moltable.molNum();
        
            if (this->mols[0].contains(molnum))
            {
                //calculate the forces on this molecule caused by group0
                int imol = this->mols[0].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols0_array[imol];
            
                for (int j=0; j<nmols1; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols1_array[j];
                
                    Potential::calculateField(mol0, mol1, probe, moltable,
                                              component, this->components(),
                                              workspace, scale_field);
                }
            }
        
            if (this->mols[1].contains(molnum))
            {
                //calculate the forces on this molecule caused by group1
                int imol = this->mols[1].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols1_array[imol];
            
                for (int j=0; j<nmols0; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols0_array[j];
                
                    Potential::calculateField(mol0, mol1, probe, moltable,
                                              component, this->components(),
                                              workspace, scale_field);
                }
            }
        }
    }

    if (ngrids > 0)
    {
        GridFieldTable *gridtable_array = fieldtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridFieldTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol0 = mols0_array[j];
                    
                Potential::calculateField(mol0, probe, gridtable, component,
                                          this->components(), workspace, scale_field);
            }

            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                    
                Potential::calculateField(mol1, probe, gridtable, component,
                                          this->components(), workspace, scale_field);
            }
        }
    }
}

/** Calculate the potential of the passed probe acting at the points in 
    the passed potentialtable
    that arise from this forcefield, and add them onto the potentials present
    in the potential table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         const Probe &prob, double scale_potential)
{
    if (scale_potential == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int ngrids = potentialtable.nGrids();
    const int npotmols = potentialtable.nMolecules();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    
    typename Potential::PotentialWorkspace workspace;
    
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();

    if (npotmols > 0)
    {
        MolPotentialTable *moltable_array = potentialtable.moleculeData();
    
        for (int i=0; i<npotmols; ++i)
        {
            MolPotentialTable &moltable = moltable_array[i];
            
            MolNum molnum = moltable.molNum();
            
            if (this->mols[0].contains(molnum))
            {
                //calculate the field on this molecule caused by group0
                int imol = this->mols[0].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols0_array[imol];
                
                for (int j=0; j<nmols1; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols1_array[j];
                    
                    Potential::calculatePotential(mol0, mol1, probe, moltable,
                                                  workspace, scale_potential);
                }
            }
            
            if (this->mols[1].contains(molnum))
            {
                //calculate the fields on this molecule caused by group1
                int imol = this->mols[1].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols1_array[imol];
                
                for (int j=0; j<nmols0; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols0_array[j];
                    
                    Potential::calculatePotential(mol0, mol1, probe, moltable,
                                                  workspace, scale_potential);
                }
            }
        }
    }
    
    if (ngrids > 0)
    {
        GridPotentialTable *gridtable_array = potentialtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridPotentialTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol0 = mols0_array[j];
                    
                Potential::calculatePotential(mol0, probe, gridtable,
                                              workspace, scale_potential);
            }

            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                    
                Potential::calculatePotential(mol1, probe, gridtable,
                                              workspace, scale_potential);
            }
        }
    }
}

/** Calculate the potential on the passed probe acting on the points in 
    the passed potential table  
    caused by the component of this forcefield represented by 'component',
    adding this potential onto the existing potentials in the potential table (optionally
    multiplied by 'scale_potential' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         const Symbol &component, const Probe &prob,
                                         double scale_potential)
{
    if (scale_potential == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int npotmols = potentialtable.nMolecules();
    const int ngrids = potentialtable.nGrids();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    
    typename Potential::PotentialWorkspace workspace;
    
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    if (npotmols > 0)
    {
        MolPotentialTable *pottable_array = potentialtable.moleculeData();

        for (int i=0; i<npotmols; ++i)
        {
            MolPotentialTable &moltable = pottable_array[i];
        
            MolNum molnum = moltable.molNum();
        
            if (this->mols[0].contains(molnum))
            {
                //calculate the forces on this molecule caused by group0
                int imol = this->mols[0].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols0_array[imol];
            
                for (int j=0; j<nmols1; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols1_array[j];
                
                    Potential::calculatePotential(mol0, mol1, probe, moltable,
                                                  component, this->components(),
                                                  workspace, scale_potential);
                }
            }
        
            if (this->mols[1].contains(molnum))
            {
                //calculate the forces on this molecule caused by group1
                int imol = this->mols[1].indexOf(molnum);
                const typename Potential::Molecule &mol0 = mols1_array[imol];
            
                for (int j=0; j<nmols0; ++j)
                {
                    const typename Potential::Molecule &mol1 = mols0_array[j];
                
                    Potential::calculatePotential(mol0, mol1, probe, moltable,
                                                  component, this->components(),
                                                  workspace, scale_potential);
                }
            }
        }
    }

    if (ngrids > 0)
    {
        GridPotentialTable *gridtable_array = potentialtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridPotentialTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol0 = mols0_array[j];
                    
                Potential::calculatePotential(mol0, probe, gridtable, component,
                                              this->components(), workspace, 
                                              scale_potential);
            }

            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                    
                Potential::calculatePotential(mol1, probe, gridtable, component,
                                              this->components(), workspace, 
                                              scale_potential);
            }
        }
    }
}

/** Calculate the fields acting at the points in the passed fieldtable
    that arise from this forcefield, and add them onto the fields present
    in the field table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::field(FieldTable &fieldtable, double scale_field)
{
    typename Potential::Probe default_probe;
    Inter2B2G3DFF<Potential>::field(fieldtable, default_probe, scale_field);
}

/** Calculate the field acting on the points in the passed forcetable  
    caused by the component of this forcefield represented by 'component',
    adding this field onto the existing fields in the forcetable (optionally
    multiplied by 'scale_field' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                     double scale_field)
{
    typename Potential::Probe default_probe;
    Inter2B2G3DFF<Potential>::field(fieldtable, component,
                                    default_probe, scale_field);
}

/** Calculate the potential acting at the points in the passed potentialtable
    that arise from this forcefield, and add them onto the potentials present
    in the potential table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         double scale_potential)
{
    typename Potential::Probe default_probe;
    Inter2B2G3DFF<Potential>::potential(potentialtable, default_probe, scale_potential);
}

/** Calculate the potential acting on the points in the passed potential table  
    caused by the component of this forcefield represented by 'component',
    adding this potential onto the existing potentials in the potential table (optionally
    multiplied by 'scale_potential' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         const Symbol &component, double scale_potential)
{
    typename Potential::Probe default_probe;
    Inter2B2G3DFF<Potential>::potential( potentialtable, component, 
                                         default_probe, scale_potential );
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireFF

SIRE_END_HEADER

#endif
