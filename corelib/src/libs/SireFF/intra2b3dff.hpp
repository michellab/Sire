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

#ifndef SIREFF_INTRA2B3DFF_HPP
#define SIREFF_INTRA2B3DFF_HPP

#include "ff3d.h"
#include "intra2bff.hpp"

SIRE_BEGIN_HEADER

namespace SireFF
{

/** This class provides an intramolecular non-bonded potential
    that can work with an two-body, three-dimensional potential
    (provided by the template class 'Potential'). This is a
    3D specific forcefield, so provides functions that can 
    calculate the 3D forces acting on the atoms.
    
    @author Christopher Woods
*/
template<class Potential>
class SIREFF_EXPORT Intra2B3DFF 
                  : public SireBase::ConcreteProperty< Intra2B3DFF<Potential>,
                                                       Intra2BFF<Potential> >, 
                    public FF3D
{
public:
    typedef typename Potential::Components Components;

    Intra2B3DFF();
    Intra2B3DFF(const QString &name);
    
    Intra2B3DFF(const Intra2B3DFF<Potential> &other);
    
    ~Intra2B3DFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Intra2B3DFF<Potential>& operator=(const Intra2B3DFF<Potential> &other);
    
    bool operator==(const Intra2B3DFF<Potential> &other) const;
    bool operator!=(const Intra2B3DFF<Potential> &other) const;
    
    Intra2B3DFF<Potential>* clone() const;
    
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
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor (without giving the forcefield a name!) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B3DFF<Potential>::Intra2B3DFF()
  : SireBase::ConcreteProperty< Intra2B3DFF<Potential>,Intra2BFF<Potential> >(), 
    FF3D()
{}

/** Construct, giving this forcefield a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B3DFF<Potential>::Intra2B3DFF(const QString &name)
  : SireBase::ConcreteProperty< Intra2B3DFF<Potential>,Intra2BFF<Potential> >(name), 
    FF3D()
{}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B3DFF<Potential>::Intra2B3DFF(const Intra2B3DFF<Potential> &other)
  : SireBase::ConcreteProperty< Intra2B3DFF<Potential>,Intra2BFF<Potential> >(other), 
    FF3D(other)
{}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B3DFF<Potential>::~Intra2B3DFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B3DFF<Potential>& 
Intra2B3DFF<Potential>::operator=(const Intra2B3DFF<Potential> &other)
{
    Intra2BFF<Potential>::operator=(other);
    FF3D::operator=(other);
    
    return *this;
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B3DFF<Potential>::operator==(const Intra2B3DFF<Potential> &other) const
{
    return Intra2BFF<Potential>::operator==(other);
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B3DFF<Potential>::operator!=(const Intra2B3DFF<Potential> &other) const
{
    return Intra2BFF<Potential>::operator!=(other);
}
    
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2B3DFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Intra2B3DFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2B3DFF<Potential>::what() const
{
    return Intra2B3DFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B3DFF<Potential>* Intra2B3DFF<Potential>::clone() const
{
    return new Intra2B3DFF<Potential>(*this);
}

/** Calculate the energies of the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the energies present
    in the energy table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "Intra2B3DFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the energies acting on the molecules in the passed energytable  
    caused by the component of this forcefield represented by 'symbol',
    adding this energy onto the existing energies in the forcetable (optionally
    multiplied by 'scale_energy' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::energy(EnergyTable &energytable, const Symbol &symbol, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "Intra2B3DFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the forces acting on the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the forces present
    in the force table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::force(ForceTable &forcetable, double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    
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
        const typename Potential::Molecule &mol = mols_array[imol];
            
        //calculate the intramolecular forces acting on this molecule
        Potential::calculateForce(mol, moltable, workspace, scale_force);
    }
}

/** Calculate the force acting on the molecules in the passed forcetable  
    caused by the component of this forcefield represented by 'symbol',
    adding this force onto the existing forces in the forcetable (optionally
    multiplied by 'scale_force' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::force(ForceTable &forcetable, const Symbol &symbol,
                                   double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    
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
        const typename Potential::Molecule &mol = mols_array[imol];
            
        //calculate the intramolecular force acting on this molecule
        Potential::calculateForce(mol, moltable, symbol,
                                  this->components(), workspace, scale_force);
    }
}

/** Calculate the fields acting at the points in the passed fieldtable
    that arise from this forcefield, and add them onto the fields present
    in the field table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::field(FieldTable &fieldtable, const Probe &prob, 
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
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<nfieldmols; ++i)
    {
        MolFieldTable &moltable = fieldtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the field
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol = mols_array[imol];
            
        //calculate the intramolecular fields acting on this molecule
        Potential::calculateField(mol, probe, moltable, workspace, scale_field);
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
void Intra2B3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
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
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<nfieldmols; ++i)
    {
        MolFieldTable &moltable = fieldtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the field
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol = mols_array[imol];
            
        //calculate the intramolecular fields acting on this molecule
        Potential::calculateField(mol, probe, moltable, 
                                  component, this->components(),
                                  workspace, scale_field);
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
void Intra2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
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
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<npotentialmols; ++i)
    {
        MolPotentialTable &moltable = potentialtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the potential
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol = mols_array[imol];
            
        //calculate the intramolecular potential acting on this molecule
        Potential::calculatePotential(mol, probe, moltable, workspace, scale_potential);
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
void Intra2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
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
    const ChunkedVector<typename Potential::Molecule> mols_array 
                            = this->mols.moleculesByIndex();
    
    for (int i=0; i<npotentialmols; ++i)
    {
        MolPotentialTable &moltable = potentialtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (not this->mols.contains(molnum))
            //we don't contain this molecule, so no point
            //calculating the potential
            continue;
            
        //get the copy of this molecule from this forcefield
        int imol = this->mols.indexOf(molnum);
        const typename Potential::Molecule &mol = mols_array[imol];
            
        //calculate the intramolecular potential acting on this molecule
        Potential::calculatePotential(mol, probe, moltable, 
                                      component, this->components(),
                                      workspace, scale_potential);
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
void Intra2B3DFF<Potential>::field(FieldTable &fieldtable, double scale_field)
{
    typename Potential::Probe default_probe;
    Intra2B3DFF<Potential>::field(fieldtable, default_probe, scale_field);
}

/** Calculate the field acting on the points in the passed forcetable  
    caused by the component of this forcefield represented by 'component',
    adding this field onto the existing fields in the forcetable (optionally
    multiplied by 'scale_field' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                   double scale_field)
{
    typename Potential::Probe default_probe;
    Intra2B3DFF<Potential>::field(fieldtable, component, default_probe, scale_field);
}

/** Calculate the potential acting at the points in the passed potentialtable
    that arise from this forcefield, and add them onto the potentials present
    in the potential table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                       double scale_potential)
{
    typename Potential::Probe default_probe;
    Intra2B3DFF<Potential>::potential(potentialtable, default_probe, scale_potential);
}

/** Calculate the potential acting on the points in the passed potential table  
    caused by the component of this forcefield represented by 'component',
    adding this potential onto the existing potentials in the potential table (optionally
    multiplied by 'scale_potential' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                       const Symbol &component, double scale_potential)
{
    typename Potential::Probe default_probe;
    Intra2B3DFF<Potential>::potential(potentialtable, component,
                                      default_probe, scale_potential);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireFF

SIRE_END_HEADER

#endif
