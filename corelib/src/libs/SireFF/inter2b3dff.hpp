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

#ifndef SIREFF_INTER2B3DFF_HPP
#define SIREFF_INTER2B3DFF_HPP

#include "ff3d.h"
#include "inter2bff.hpp"

#include "SireBase/countflops.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireFF
{

/** This class provides an intermolecular non-bonded potential
    that can work with an two-body, three-dimensional potential
    (provided by the template class 'Potential'). This is a
    3D specific forcefield, so provides functions that can 
    calculate the 3D forces acting on the atoms.
    
    @author Christopher Woods
*/
template<class Potential>
class Inter2B3DFF 
              : public SireBase::ConcreteProperty< Inter2B3DFF<Potential>,
                                                   Inter2BFF<Potential> >, 
                public FF3D
{
public:
    Inter2B3DFF();
    Inter2B3DFF(const QString &name);
    
    Inter2B3DFF(const Inter2B3DFF<Potential> &other);
    
    ~Inter2B3DFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Inter2B3DFF<Potential>& operator=(const Inter2B3DFF<Potential> &other);
    
    bool operator==(const Inter2B3DFF<Potential> &other) const;
    bool operator!=(const Inter2B3DFF<Potential> &other) const;
    
    Inter2B3DFF<Potential>* clone() const;

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

    void packCoordinates();

protected:
    typedef typename Inter2BFF<Potential>::Energy Energy;
    typedef typename Inter2BFF<Potential>::EnergyWorkspace EnergyWorkspace;

    typedef typename Inter2BFF<Potential>::Molecules Molecules;
    typedef typename Inter2BFF<Potential>::Molecule Molecule;

    typedef typename Inter2BFF<Potential>::ChangedMolecule ChangedMolecule;

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
Inter2B3DFF<Potential>::Inter2B3DFF()
  : SireBase::ConcreteProperty< Inter2B3DFF<Potential>,Inter2BFF<Potential> >(), 
    FF3D()
{}

/** Construct, giving this forcefield a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B3DFF<Potential>::Inter2B3DFF(const QString &name)
  : SireBase::ConcreteProperty< Inter2B3DFF<Potential>,Inter2BFF<Potential> >(name), 
    FF3D()
{}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B3DFF<Potential>::Inter2B3DFF(const Inter2B3DFF<Potential> &other)
  : SireBase::ConcreteProperty< Inter2B3DFF<Potential>,Inter2BFF<Potential> >(other), 
    FF3D(other)
{}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B3DFF<Potential>::~Inter2B3DFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B3DFF<Potential>& 
Inter2B3DFF<Potential>::operator=(const Inter2B3DFF<Potential> &other)
{
    Inter2BFF<Potential>::operator=(other);
    FF3D::operator=(other);
    
    return *this;
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2B3DFF<Potential>::operator==(const Inter2B3DFF<Potential> &other) const
{
    return Inter2BFF<Potential>::operator==(other);
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2B3DFF<Potential>::operator!=(const Inter2B3DFF<Potential> &other) const
{
    return Inter2BFF<Potential>::operator!=(other);
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Inter2B3DFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Inter2B3DFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Inter2B3DFF<Potential>::what() const
{
    return Inter2B3DFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2B3DFF<Potential>* Inter2B3DFF<Potential>::clone() const
{
    return new Inter2B3DFF<Potential>(*this);
}

/** Pack the coordinates of the molecules so that they all lie 
    contiguously in memory */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::packCoordinates()
{
    this->mols.packCoordinates();
}

/** Recalculate the energy of the current state of this forcefield. This
    will recalculate the energy using the quickest possible route, e.g.
    if will only recalculate the energies of molecules that have changed
    since the last evaluation */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::recalculateEnergy()
{
    int nmols = this->mols.count();
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();

    if (this->changed_mols.count() == nmols)
        //all of the molecules have changed!
        this->changed_mols.clear();

    //tell the potential that we are starting an evaluation
    Potential::startEvaluation();

    try
    {

    if (this->changed_mols.isEmpty())
    {
        Energy total_nrg;

        {
            //we are not recording changes, so we have to assume that
            //everything has changed. Recalculate the total energy from scratch
            EnergyWorkspace workspace;
            Energy my_total_nrg;

            const ChunkedVector<typename Potential::Molecule> &my_mols_array = mols_array;
            const ChunkedVector<SireVol::AABox> &aaboxes_array 
                                                = this->mols.aaBoxesByIndex();
            const int my_nmols = nmols;

            const SireVol::Space &spce = this->space();
            const double cutoff = this->switchingFunction().cutoffDistance();

            //loop over all pairs of molecules
            for (int i=0; i<my_nmols-1; ++i)
            {
                const typename Potential::Molecule &mol0 = my_mols_array.at(i);
                const SireVol::AABox &aabox0 = aaboxes_array.at(i);
        
                for (int j=i+1; j<my_nmols; ++j)
                {
                    if (not spce.beyond(cutoff, aabox0, aaboxes_array.at(j)))
                    {
                        const typename Potential::Molecule &mol1 = my_mols_array.at(j);
                        Potential::calculateEnergy(mol0, mol1, my_total_nrg, workspace);
                    }
                }
            }
              
            {
                total_nrg += my_total_nrg;
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

        {
            EnergyWorkspace old_workspace;
            EnergyWorkspace new_workspace;

            Energy my_old_nrg;
            Energy my_new_nrg;

            const QHash<MolNum,ChangedMolecule> my_changed_mols = this->changed_mols;
            const ChunkedVector<typename Potential::Molecule> &my_mols_array = mols_array;
            const ChunkedVector<SireVol::AABox> aaboxes_array 
                                                    = this->mols.aaBoxesByIndex();
            const int my_nmols = nmols;

            const SireVol::Space &spce = this->space();
            const double cutoff = this->switchingFunction().cutoffDistance();

            if (my_changed_mols.count() == 1)
            {
                const typename Potential::ChangedMolecule &changed_mol = 
                                        *(my_changed_mols.constBegin());
                         
                int idx = this->mols.indexOf(changed_mol.number());

                const SireVol::AABox &new_box = changed_mol.newMolecule().aaBox();
                const SireVol::AABox &old_box = changed_mol.oldMolecule().aaBox();

                for (int i=0; i<my_nmols; ++i)
                {
                    if (idx != i)
                    {
                        const typename Potential::Molecule &mol = my_mols_array.at(i);
                        const SireVol::AABox &aabox = aaboxes_array.at(i);
                    
                        if (not spce.beyond(cutoff, old_box, aabox))
                        {
                            Potential::calculateEnergy(mol, changed_mol.oldParts(),
                                                       my_old_nrg, old_workspace);
                        }
                        
                        if (not spce.beyond(cutoff, new_box, aabox))
                        {
                            Potential::calculateEnergy(mol, changed_mol.newParts(),
                                                       my_new_nrg, new_workspace);
                        }
                    }
                }
            }
            else
            {
                for (int i=0; i<my_nmols; ++i)
                {
                    const typename Potential::Molecule &mol0 = my_mols_array.at(i);
                
                    typename QHash<MolNum,ChangedMolecule>::const_iterator it
                                             = my_changed_mols.constFind(mol0.number());
                    
                    //get the bounding box for this molecule
                    const SireVol::AABox &aabox0 = aaboxes_array.at(i);
                                                                                                                 
                    if (it == my_changed_mols.constEnd())
                    {
                        //this molecule has not changed - just calculate its
                        //energy with all of the changed molecules - as this molecule
                        //hasn't changed, we only need to calculate the change
                        //in energy between this molecule and the changed parts
                        //of the changed molecules
                        for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                          it2 = my_changed_mols.constBegin();
                             it2 != this->changed_mols.constEnd();
                             ++it2)
                        {
                            if (not spce.beyond(cutoff, aabox0, it2->oldMolecule().aaBox()))
                            {
                                Potential::calculateEnergy(mol0, it2->oldParts(),
                                                           my_old_nrg, old_workspace);
                            }
                        
                            if (not spce.beyond(cutoff, aabox0, it2->newMolecule().aaBox()))
                            {
                                Potential::calculateEnergy(mol0, it2->newParts(),
                                                           my_new_nrg, new_workspace);
                            }
                        }
                    }
                    else if (this->changed_mols.count() > 1)
                    {
                        //this molecule has changed - calculate its energy with all
                        //of the changed molecules that lie after it in the changed_mols
                        //hash (thus ensuring we don't double-count)
                        typename QHash<MolNum,ChangedMolecule>::const_iterator it2 = it;
                    
                        bool this_changed_all = it->changedAll();
                    
                        if (this_changed_all)
                        {
                            //all of this molecule has changed - so we need to 
                            //calculate the energy of this molecule with *all* of
                            //the parts of the other changed molecules
                            for (++it2; it2 != my_changed_mols.constEnd(); ++it2)
                            {
                                Potential::calculateEnergy(it->oldMolecule(),
                                                           it2->oldMolecule(),
                                                           my_old_nrg, old_workspace);
                                                       
                                Potential::calculateEnergy(it->newMolecule(),
                                                           it2->newMolecule(),
                                                           my_new_nrg, new_workspace);
                            }
                        }
                        else
                        {
                            for (++it2; it2 != my_changed_mols.constEnd(); ++it2)
                            {
                                if (it2->changedAll())
                                {
                                    //all of the other molecule has changed - we need
                                    //to calculate the energy of the whole molecules
                                    //interaction
                                    Potential::calculateEnergy(it->oldMolecule(),
                                                               it2->oldMolecule(),
                                                               my_old_nrg, old_workspace);
                                
                                    Potential::calculateEnergy(it->newMolecule(),
                                                               it2->newMolecule(),
                                                               my_new_nrg, new_workspace);
                                }
                                else
                                {
                                    //both a part of this molecule and a part of the 
                                    //other molecule have changed
                               
                                    //the change in energy associated with changing 
                                    //the first molecule...
                                    Potential::calculateEnergy(it->oldParts(),
                                                               it2->oldMolecule(),
                                                               my_old_nrg, old_workspace);
                                                           
                                    Potential::calculateEnergy(it->newParts(),
                                                               it2->newMolecule(),
                                                               my_new_nrg, new_workspace);
                                                           
                                    //now the change in energy associated with changing
                                    //the second molecule...
                                    Potential::calculateEnergy(it2->oldParts(),
                                                               it->oldMolecule(),
                                                               my_old_nrg, old_workspace);
                                                         
                                    Potential::calculateEnergy(it2->newParts(),
                                                               it->newMolecule(),
                                                               my_new_nrg, new_workspace);
                                                          
                                    //now remove double counted changed in mol1 with
                                    //change in mol2
                                    Potential::calculateEnergy(it->oldParts(),
                                                               it2->oldParts(),
                                                               my_old_nrg, old_workspace, -1);
                                                          
                                    Potential::calculateEnergy(it->newParts(),
                                                               it2->newParts(),
                                                               my_new_nrg, new_workspace, -1);
                                }
                            }
                        }
                    }
                }
            }

            {
                old_nrg += my_old_nrg;
                new_nrg += my_new_nrg;
            }

        } // end of parallel block

        if ( not this->removed_mols.isEmpty() )
        {
            //finally, loop over all of the molecules that have been removed - the energy
            //of non-changed molecules with removed molecules has already been calculated,
            //as has the energy of moved molecules that are before the removed molecules
            //in the moved list. We only now have to calculate the energy of the removed
            //molecules with all of the molecules that lie above us in the moved list
            
            EnergyWorkspace workspace;
            
            for (typename QSet<MolNum>::const_iterator 
                                                it = this->removed_mols.constBegin();
                 it != this->removed_mols.constEnd();
                 ++it)
            {
                typename QSet<MolNum>::const_iterator it2 = it;
            
                const ChangedMolecule &mol0 = *(this->changed_mols.constFind(*it));
            
                for (++it2; it2 != this->removed_mols.constEnd(); ++it2)
                {
                    const ChangedMolecule &mol1 = *(this->changed_mols.constFind(*it2));
            
                    Potential::calculateEnergy(mol0.oldMolecule(),
                                               mol1.oldMolecule(),
                                               old_nrg, workspace);
                                           
                    //molecule has been removed, so no new energy
                }
            }
        }
         
        //change the energy
        this->components().changeEnergy(*this, new_nrg - old_nrg);
        
        //clear the changed molecules
        this->changed_mols.clear();
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
SireUnits::Dimension::MolarEnergy Inter2B3DFF<Potential>::energy()
{
    return Inter2BFF<Potential>::energy();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
SireUnits::Dimension::MolarEnergy Inter2B3DFF<Potential>::energy(const Symbol &component)
{
    return Inter2BFF<Potential>::energy(component);
}

/** Calculate the energies of the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the energies present
    in the energy table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::energy(EnergyTable &energytable, double scale_energy)
{

    if (scale_energy == 0)
        return;

    int nenergymols = energytable.count();
    int nmols = this->mols.count();
    
    typename Potential::EnergyWorkspace workspace;
    
    MolEnergyTable *energytable_array = energytable.data();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = this->mols.moleculesByIndex();
    
    typename Potential::Energy energy;

    for (int i=0; i<nenergymols; ++i)
    {
        MolEnergyTable &moltable = energytable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the force acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculateEnergy(mol0, mol1, moltable, workspace, scale_energy);
        }
    }
}

/** Calculate the energies acting on the molecules in the passed energytable  
    caused by the component of this forcefield represented by 'symbol',
    adding this energy onto the existing energies in the forcetable (optionally
    multiplied by 'scale_energy' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::energy(EnergyTable &energytable, const Symbol &symbol, double scale_energy)
{
    if (scale_energy == 0)
        return;

    int nenergymols = energytable.count();
    int nmols = this->mols.count();
    
    typename Potential::EnergyWorkspace workspace;
    
    MolEnergyTable *energytable_array = energytable.data();
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();
    
    typename Potential::Energy energy;

    for (int i=0; i<nenergymols; ++i)
    {
        MolEnergyTable &moltable = energytable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the force acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];

            Potential::calculateEnergy(mol0, mol1, moltable, workspace, scale_energy);

        }
    }

 }

/** Calculate the forces acting on the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the forces present
    in the force table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::force(ForceTable &forcetable, double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    int nmols = this->mols.count();
    
    typename Potential::ForceWorkspace workspace;
    
    MolForceTable *forcetable_array = forcetable.data();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the force acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculateForce(mol0, mol1, moltable, workspace, scale_force);
        }
    }
}

/** Calculate the force acting on the molecules in the passed forcetable  
    caused by the component of this forcefield represented by 'symbol',
    adding this force onto the existing forces in the forcetable (optionally
    multiplied by 'scale_force' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::force(ForceTable &forcetable, const Symbol &symbol,
                                   double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    int nmols = this->mols.count();
    
    typename Potential::ForceWorkspace workspace;
    
    MolForceTable *forcetable_array = forcetable.data();
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the force acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculateForce(mol0, mol1, moltable, symbol,
                                      this->components(), workspace, scale_force);
        }
    }
}

/** Calculate the fields acting at the points in the passed fieldtable
    that arise from this forcefield, and add them onto the fields present
    in the field table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::field(FieldTable &fieldtable, const Probe &prob, 
                                   double scale_field)
{
    if (scale_field == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int ngrids = fieldtable.nGrids();
    const int nfieldmols = fieldtable.nMolecules();
    const int nmols = this->mols.count();
    
    typename Potential::FieldWorkspace workspace;
    
    MolFieldTable *fieldtable_array = fieldtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<nfieldmols; ++i)
    {
        MolFieldTable &moltable = fieldtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the field acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculateField(mol0, mol1, probe, moltable, 
                                      workspace, scale_field);
        }
    }

    if (ngrids > 0)
    {
        GridFieldTable *gridtable_array = fieldtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridFieldTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols; ++j)
            {
                const typename Potential::Molecule &mol = mols_array[j];
                    
                Potential::calculateField(mol, probe, gridtable,
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
void Inter2B3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                   const Probe &prob, double scale_field)
{
    if (scale_field == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int ngrids = fieldtable.nGrids();
    const int nfieldmols = fieldtable.nMolecules();
    const int nmols = this->mols.count();
    
    typename Potential::FieldWorkspace workspace;
    
    MolFieldTable *fieldtable_array = fieldtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<nfieldmols; ++i)
    {
        MolFieldTable &moltable = fieldtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the field acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculateField(mol0, mol1, probe, moltable,
                                      component, this->components(), 
                                      workspace, scale_field);
        }
    }

    if (ngrids > 0)
    {
        GridFieldTable *gridtable_array = fieldtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridFieldTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols; ++j)
            {
                const typename Potential::Molecule &mol = mols_array[j];
                    
                Potential::calculateField(mol, probe, gridtable,
                                          component, this->components(),
                                          workspace, scale_field);
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
void Inter2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                       const Probe &prob, double scale_potential)
{
    if (scale_potential == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int ngrids = potentialtable.nGrids();
    const int npotentialmols = potentialtable.nMolecules();
    const int nmols = this->mols.count();
    
    typename Potential::PotentialWorkspace workspace;
    
    MolPotentialTable *potentialtable_array = potentialtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<npotentialmols; ++i)
    {
        MolPotentialTable &moltable = potentialtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the field acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculatePotential(mol0, mol1, probe, moltable, 
                                          workspace, scale_potential);
        }
    }

    if (ngrids > 0)
    {
        GridPotentialTable *gridtable_array = potentialtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridPotentialTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols; ++j)
            {
                const typename Potential::Molecule &mol = mols_array[j];
                    
                Potential::calculatePotential(mol, probe, gridtable,
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
void Inter2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                       const Symbol &component,
                                       const Probe &prob, double scale_potential)
{
    if (scale_potential == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int ngrids = potentialtable.nGrids();
    const int npotentialmols = potentialtable.nMolecules();
    const int nmols = this->mols.count();
    
    typename Potential::PotentialWorkspace workspace;
    
    MolPotentialTable *potentialtable_array = potentialtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<npotentialmols; ++i)
    {
        MolPotentialTable &moltable = potentialtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the force
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol0 = mols_array[imol];
            
        //calculate the field acting on this molecule caused by all of the 
        //other molecules in this forcefield
        for (int j=0; j<nmols; ++j)
        {
            if (j == imol)
                continue;
                
            const typename Potential::Molecule &mol1 = mols_array[j];
            
            Potential::calculatePotential(mol0, mol1, probe, moltable,
                                          component, this->components(),
                                          workspace, scale_potential);
        }
    }

    if (ngrids > 0)
    {
        GridPotentialTable *gridtable_array = potentialtable.gridData();
        
        for (int i=0; i<ngrids; ++i)
        {
            GridPotentialTable &gridtable = gridtable_array[i];
            
            for (int j=0; j<nmols; ++j)
            {
                const typename Potential::Molecule &mol = mols_array[j];
                    
                Potential::calculatePotential(mol, probe, gridtable,
                                              component, this->components(),
                                              workspace, scale_potential);
            }
        }
    }
}

/** Calculate the fields acting at the points in the passed fieldtable
    that arise from this forcefield, and add them onto the fields present
    in the field table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::field(FieldTable &fieldtable, double scale_field)
{
    typename Potential::Probe default_probe;
    Inter2B3DFF<Potential>::field(fieldtable, default_probe, scale_field);
}

/** Calculate the field acting on the points in the passed forcetable  
    caused by the component of this forcefield represented by 'component',
    adding this field onto the existing fields in the forcetable (optionally
    multiplied by 'scale_field' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                   double scale_field)
{
    typename Potential::Probe default_probe;
    Inter2B3DFF<Potential>::field(fieldtable, component, default_probe, scale_field);
}

/** Calculate the potential acting at the points in the passed potentialtable
    that arise from this forcefield, and add them onto the potentials present
    in the potential table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                       double scale_potential)
{
    typename Potential::Probe default_probe;
    Inter2B3DFF<Potential>::potential(potentialtable, default_probe, scale_potential);
}

/** Calculate the potential acting on the points in the passed potential table  
    caused by the component of this forcefield represented by 'component',
    adding this potential onto the existing potentials in the potential table (optionally
    multiplied by 'scale_potential' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                       const Symbol &component,
                                       double scale_potential)
{
    typename Potential::Probe default_probe;
    Inter2B3DFF<Potential>::potential(potentialtable, component, 
                                      default_probe, scale_potential);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireFF

SIRE_END_HEADER

#endif
