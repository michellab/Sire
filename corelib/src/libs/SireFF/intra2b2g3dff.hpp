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

#ifndef SIREFF_INTRA2B2G3DFF_HPP
#define SIREFF_INTRA2B2G3DFF_HPP

#include "ff3d.h"
#include "intra2b2gff.hpp"

SIRE_BEGIN_HEADER

namespace SireFF
{

/** This class provides an intramolecular non-bonded potential
    that can work with an two-body, three-dimensional potential
    (provided by the template class 'Potential') that
    calculates the intramolecular potential between
    two groups of molecules. This is a 3D specific forcefield, 
    so provides functions that can calculate the 3D forces acting on the atoms.
    
    (so you could use this to calculate the forces acting on one
    part of a molecule that are caused by another part of the molecule)
    
    @author Christopher Woods
*/
template<class Potential>
class SIREFF_EXPORT Intra2B2G3DFF 
                    : public SireBase::ConcreteProperty< Intra2B2G3DFF<Potential>,
                                                         Intra2B2GFF<Potential> >, 
                      public FF3D
{
public:
    Intra2B2G3DFF();
    Intra2B2G3DFF(const QString &name);
    
    Intra2B2G3DFF(const Intra2B2G3DFF<Potential> &other);
    
    ~Intra2B2G3DFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Intra2B2G3DFF<Potential>& operator=(const Intra2B2G3DFF<Potential> &other);
    
    bool operator==(const Intra2B2G3DFF<Potential> &other) const;
    bool operator!=(const Intra2B2G3DFF<Potential> &other) const;
    
    Intra2B2G3DFF<Potential>* clone() const;

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
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor (without giving the forcefield a name!) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2G3DFF<Potential>::Intra2B2G3DFF()
  : SireBase::ConcreteProperty< Intra2B2G3DFF<Potential>,Intra2B2GFF<Potential> >(), 
    FF3D()
{}

/** Construct, giving this forcefield a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2G3DFF<Potential>::Intra2B2G3DFF(const QString &name)
  : SireBase::ConcreteProperty< Intra2B2G3DFF<Potential>,Intra2B2GFF<Potential> >(name), 
    FF3D()
{}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2G3DFF<Potential>::Intra2B2G3DFF(const Intra2B2G3DFF<Potential> &other)
  : SireBase::ConcreteProperty< Intra2B2G3DFF<Potential>,Intra2B2GFF<Potential> >(other), 
    FF3D(other)
{}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2G3DFF<Potential>::~Intra2B2G3DFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2G3DFF<Potential>& 
Intra2B2G3DFF<Potential>::operator=(const Intra2B2G3DFF<Potential> &other)
{
    Intra2B2GFF<Potential>::operator=(other);
    FF3D::operator=(other);
    
    return *this;
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2G3DFF<Potential>::operator==(const Intra2B2G3DFF<Potential> &other) const
{
    return Intra2B2GFF<Potential>::operator==(other);
}

/** Comparison operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2G3DFF<Potential>::operator!=(const Intra2B2G3DFF<Potential> &other) const
{
    return Intra2B2GFF<Potential>::operator!=(other);
}
    
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2B2G3DFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Intra2B2G3DFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2B2G3DFF<Potential>::what() const
{
    return Intra2B2G3DFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2G3DFF<Potential>* Intra2B2G3DFF<Potential>::clone() const
{
    return new Intra2B2G3DFF<Potential>(*this);
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
SireUnits::Dimension::MolarEnergy Intra2B2G3DFF<Potential>::energy()
{
    return Intra2B2GFF<Potential>::energy();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
SireUnits::Dimension::MolarEnergy Intra2B2G3DFF<Potential>::energy(const Symbol &component)
{
    return Intra2B2GFF<Potential>::energy(component);
}

/** Calculate the energies of the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the energies present
    in the energy table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "Intra2B2G3DFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the energies acting on the molecules in the passed energytable  
    caused by the component of this forcefield represented by 'symbol',
    adding this energy onto the existing energies in the forcetable (optionally
    multiplied by 'scale_energy' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::energy(EnergyTable &energytable, const Symbol &symbol, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "Intra2B2G3DFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the forces acting on the molecules in the passed forcetable
    that arise from this forcefield, and add them onto the forces present
    in the force table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::force(ForceTable &forcetable, double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    
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
        
        if (this->mols[0].contains(molnum) and this->mols[1].contains(molnum))
        {
            int imol0 = this->mols[0].indexOf(molnum);
            int imol1 = this->mols[1].indexOf(molnum);
            
            const typename Potential::Molecule &mol0 = mols0_array[imol0];
            const typename Potential::Molecule &mol1 = mols1_array[imol1];

            //calculate the forces on mols[0] caused by mols[1]
            Potential::calculateForce(mol0, mol1, moltable, 
                                      workspace, scale_force);
                                      
            //now add on the forces on mols[1] caused by mols[0]
            Potential::calculateForce(mol1, mol0, moltable,
                                      workspace, scale_force);
        }
    }
}

/** Calculate the force acting on the molecules in the passed forcetable  
    caused by the component of this forcefield represented by 'symbol',
    adding this force onto the existing forces in the forcetable (optionally
    multiplied by 'scale_force' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::force(ForceTable &forcetable, const Symbol &symbol,
                                   double scale_force)
{
    if (scale_force == 0)
        return;

    int nforcemols = forcetable.count();
    
    typename Potential::ForceWorkspace workspace;
    
    MolForceTable *forcetable_array = forcetable.data();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<nforcemols; ++i)
    {
        MolForceTable &moltable = forcetable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum) and this->mols[1].contains(molnum))
        {
            int imol0 = this->mols[0].indexOf(molnum);
            int imol1 = this->mols[1].indexOf(molnum);
            
            const typename Potential::Molecule &mol0 = mols0_array[imol0];
            const typename Potential::Molecule &mol1 = mols1_array[imol1];

            //calculate the forces on mols[0] caused by mols[1]
            Potential::calculateForce(mol0, mol1, moltable,
                                      symbol, this->components(),
                                      workspace, scale_force);
                                      
            //now add on the forces on mols[1] caused by mols[0]
            Potential::calculateForce(mol1, mol0, moltable,
                                      symbol, this->components(),
                                      workspace, scale_force);
        }
    }
}

/** Calculate the fields acting at the points in the passed fieldtable
    that arise from this forcefield, and add them onto the fields present
    in the field table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::field(FieldTable &fieldtable, const Probe &prob, 
                                     double scale_field)
{
    if (scale_field == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int nfieldmols = fieldtable.nMolecules();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    const int ngrids = fieldtable.nGrids();
    
    typename Potential::FieldWorkspace workspace;
    
    MolFieldTable *fieldtable_array = fieldtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<nfieldmols; ++i)
    {
        MolFieldTable &moltable = fieldtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum) and this->mols[1].contains(molnum))
        {
            int imol0 = this->mols[0].indexOf(molnum);
            int imol1 = this->mols[1].indexOf(molnum);
            
            const typename Potential::Molecule &mol0 = mols0_array[imol0];
            const typename Potential::Molecule &mol1 = mols1_array[imol1];

            //calculate the field on mols[0] caused by mols[1]
            Potential::calculateField(mol0, mol1, probe, moltable, 
                                      workspace, scale_field);
                                      
            //now add on the forces on mols[1] caused by mols[0]
            Potential::calculateField(mol1, mol0, probe, moltable,
                                      workspace, scale_field);
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
                const typename Potential::Molecule &mol = mols0_array[j];
                    
                Potential::calculateField(mol, probe, gridtable,
                                          workspace, scale_field);
            }
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol = mols1_array[j];
                    
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
void Intra2B2G3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                     const Probe &prob, double scale_field)
{
    if (scale_field == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int nfieldmols = fieldtable.nMolecules();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    const int ngrids = fieldtable.nGrids();
    
    typename Potential::FieldWorkspace workspace;
    
    MolFieldTable *fieldtable_array = fieldtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<nfieldmols; ++i)
    {
        MolFieldTable &moltable = fieldtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum) and this->mols[1].contains(molnum))
        {
            int imol0 = this->mols[0].indexOf(molnum);
            int imol1 = this->mols[1].indexOf(molnum);
            
            const typename Potential::Molecule &mol0 = mols0_array[imol0];
            const typename Potential::Molecule &mol1 = mols1_array[imol1];

            //calculate the field on mols[0] caused by mols[1]
            Potential::calculateField(mol0, mol1, probe, moltable, 
                                      component, this->components(),
                                      workspace, scale_field);
                                      
            //now add on the forces on mols[1] caused by mols[0]
            Potential::calculateField(mol1, mol0, probe, moltable,
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
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol = mols0_array[j];
                    
                Potential::calculateField(mol, probe, gridtable,
                                          component, this->components(),
                                          workspace, scale_field);
            }
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol = mols1_array[j];
                    
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
void Intra2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         const Probe &prob,
                                         double scale_potential)
{
    if (scale_potential == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int npotentialmols = potentialtable.nMolecules();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    const int ngrids = potentialtable.nGrids();
    
    typename Potential::PotentialWorkspace workspace;
    
    MolPotentialTable *potentialtable_array = potentialtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<npotentialmols; ++i)
    {
        MolPotentialTable &moltable = potentialtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum) and this->mols[1].contains(molnum))
        {
            int imol0 = this->mols[0].indexOf(molnum);
            int imol1 = this->mols[1].indexOf(molnum);
            
            const typename Potential::Molecule &mol0 = mols0_array[imol0];
            const typename Potential::Molecule &mol1 = mols1_array[imol1];

            //calculate the potential on mols[0] caused by mols[1]
            Potential::calculatePotential(mol0, mol1, probe, moltable, 
                                          workspace, scale_potential);
                                      
            //now add on the potential on mols[1] caused by mols[0]
            Potential::calculatePotential(mol1, mol0, probe, moltable,
                                          workspace, scale_potential);
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
                const typename Potential::Molecule &mol = mols0_array[j];
                    
                Potential::calculatePotential(mol, probe, gridtable,
                                              workspace, scale_potential);
            }
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol = mols1_array[j];
                    
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
void Intra2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         const Symbol &component,
                                         const Probe &prob, double scale_potential)
{
    if (scale_potential == 0)
        return;

    const typename Potential::Probe probe = prob.convertTo<typename Potential::Probe>();

    const int npotentialmols = potentialtable.nMolecules();
    const int nmols0 = this->mols[0].count();
    const int nmols1 = this->mols[1].count();
    const int ngrids = potentialtable.nGrids();
    
    typename Potential::PotentialWorkspace workspace;
    
    MolPotentialTable *potentialtable_array = potentialtable.moleculeData();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = this->mols[0].moleculesByIndex();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = this->mols[1].moleculesByIndex();
    
    for (int i=0; i<npotentialmols; ++i)
    {
        MolPotentialTable &moltable = potentialtable_array[i];
        
        MolNum molnum = moltable.molNum();
        
        if (this->mols[0].contains(molnum) and this->mols[1].contains(molnum))
        {
            int imol0 = this->mols[0].indexOf(molnum);
            int imol1 = this->mols[1].indexOf(molnum);
            
            const typename Potential::Molecule &mol0 = mols0_array[imol0];
            const typename Potential::Molecule &mol1 = mols1_array[imol1];

            //calculate the potential on mols[0] caused by mols[1]
            Potential::calculatePotential(mol0, mol1, probe, moltable, 
                                          component, this->components(),
                                          workspace, scale_potential);
                                      
            //now add on the potential on mols[1] caused by mols[0]
            Potential::calculatePotential(mol1, mol0, probe, moltable,
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
            
            for (int j=0; j<nmols0; ++j)
            {
                const typename Potential::Molecule &mol = mols0_array[j];
                    
                Potential::calculatePotential(mol, probe, gridtable,
                                              component, this->components(),
                                              workspace, scale_potential);
            }
            
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol = mols1_array[j];
                    
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
void Intra2B2G3DFF<Potential>::field(FieldTable &fieldtable, double scale_field)
{
    typename Potential::Probe default_probe;
    Intra2B2G3DFF<Potential>::field(fieldtable, default_probe, scale_field);
}

/** Calculate the field acting on the points in the passed forcetable  
    caused by the component of this forcefield represented by 'component',
    adding this field onto the existing fields in the forcetable (optionally
    multiplied by 'scale_field' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::field(FieldTable &fieldtable, const Symbol &component,
                                     double scale_field)
{
    typename Potential::Probe default_probe;
    Intra2B2G3DFF<Potential>::field(fieldtable, component, default_probe, scale_field);
}

/** Calculate the potential acting at the points in the passed potentialtable
    that arise from this forcefield, and add them onto the potentials present
    in the potential table, multiplied by the passed (optional) scaling factor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         double scale_potential)
{
    typename Potential::Probe default_probe;
    Intra2B2G3DFF<Potential>::potential(potentialtable, default_probe, scale_potential);
}

/** Calculate the potential acting on the points in the passed potential table  
    caused by the component of this forcefield represented by 'component',
    adding this potential onto the existing potentials in the potential table (optionally
    multiplied by 'scale_potential' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2G3DFF<Potential>::potential(PotentialTable &potentialtable, 
                                         const Symbol &component, double scale_potential)
{
    typename Potential::Probe default_probe;
    Intra2B2G3DFF<Potential>::potential(potentialtable, component,
                                        default_probe, scale_potential);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireFF

SIRE_END_HEADER

#endif
