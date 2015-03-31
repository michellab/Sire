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

#ifndef SIREFF_INTRA2B2GFF_HPP
#define SIREFF_INTRA2B2GFF_HPP

#include "ff.h"
#include "g2ff.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
template<class Potential>
class Intra2B2GFF;
}

template<class Potential>
QDataStream& operator<<(QDataStream&, const SireFF::Intra2B2GFF<Potential>&);

template<class Potential>
QDataStream& operator>>(QDataStream&, SireFF::Intra2B2GFF<Potential>&);

namespace SireFF
{

namespace detail
{
void throwIntra2B2GFFIncompatibeScaleFactorsError(const QString &molname,
                     MolNum molnum, quint64 v0, quint64 v1);
}

/** This class provides an intramolecular non-bonded forcefield
    that can work with any two-body potential (provided via the 
    template type 'Potential') that calculates the intramolecular
    energy between two groups of molecules (this doesn't calculate
    the intramolecular energy within either group, just between
    the two groups). Use this to calculate the intramolecular energy
    between two parts of a molecule 
    
    @author Christopher Woods
*/
template<class Potential>
class SIREFF_EXPORT Intra2B2GFF 
                  : public SireBase::ConcreteProperty<Intra2B2GFF<Potential>, G2FF>, 
                    public Potential
{

friend QDataStream& ::operator<<<>(QDataStream&, const Intra2B2GFF<Potential>&);
friend QDataStream& ::operator>><>(QDataStream&, Intra2B2GFF<Potential>&);

public:
    typedef typename Potential::Components Components;

    Intra2B2GFF();
    Intra2B2GFF(const QString &name);
    
    Intra2B2GFF(const Intra2B2GFF<Potential> &other);
    
    ~Intra2B2GFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Intra2B2GFF<Potential>& operator=(const Intra2B2GFF<Potential> &other);
    
    bool operator==(const Intra2B2GFF<Potential> &other) const;
    bool operator!=(const Intra2B2GFF<Potential> &other) const;

    Intra2B2GFF<Potential>* clone() const;

    const Components& components() const;

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    void mustNowRecalculateFromScratch();    

protected:

    typedef typename Potential::Energy Energy;
    typedef typename Potential::EnergyWorkspace EnergyWorkspace;

    typedef typename Potential::Molecules Molecules;
    typedef typename Potential::Molecule Molecule;

    typedef typename Potential::ChangedMolecule ChangedMolecule;

    void assertSameIntraScaleFactors(MolNum molnum) const;

    bool recordingChanges() const;

    void recordChange(quint32 groupid, const ChangedMolecule &change);

    const FFComponent& _pvt_components() const;

    void recalculateEnergy();

    void changedPotential();
    
    void _pvt_added(quint32 groupid, const PartialMolecule &mol, 
                    const PropertyMap &map);

    void _pvt_removed(quint32 groupid, const PartialMolecule &mol);

    void _pvt_changed(quint32 groupid, const SireMol::Molecule &molecule, bool auto_commit);
    void _pvt_changed(quint32 groupid, const QList<SireMol::Molecule> &molecules, bool auto_commit);
    
    void _pvt_removedAll(quint32 groupid);
        
    bool _pvt_wouldChangeProperties(quint32 groupid, MolNum molnum, 
                                    const PropertyMap &map) const;

    void _pvt_updateName();

    /** The parameterised version of the molecules in this forcefield
        for the two groups */
    Molecules mols[2];

    /** The list of molecules that have changed since the last evaluation.
        While ffmols only contains the newest version of the molecule,
        this list contains both the newest version, and the version of the
        molecule at the last energy evaluation. */
    QHash<MolNum,ChangedMolecule> changed_mols[2];

    /** The energy components available for this forcefield */
    Components ffcomponents;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2GFF<Potential>::Intra2B2GFF() 
                       : SireBase::ConcreteProperty<Intra2B2GFF<Potential>,G2FF>(), 
                         Potential()
{
    this->_pvt_updateName();
}

/** Construct this forcefield, providing it with a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2GFF<Potential>::Intra2B2GFF(const QString &ffname)
                       : SireBase::ConcreteProperty<Intra2B2GFF<Potential>,G2FF>(), 
                         Potential()
{
    FF::setName(ffname);
}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2GFF<Potential>::Intra2B2GFF(const Intra2B2GFF<Potential> &other)
                     : SireBase::ConcreteProperty<Intra2B2GFF<Potential>,G2FF>(other), 
                       Potential(other),
                       ffcomponents(other.ffcomponents)
{
    for (int i=0; i<2; ++i)
    {
        mols[i] = other.mols[i];
        changed_mols[i] = other.changed_mols[i];
    }
}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2GFF<Potential>::~Intra2B2GFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2GFF<Potential>& 
Intra2B2GFF<Potential>::operator=(const Intra2B2GFF<Potential> &other)
{
    if (this != &other)
    {
        G2FF::operator=(other);
        Potential::operator=(other);
        
        for (int i=0; i<2; ++i)
        {
            mols[i] = other.mols[i];
            changed_mols[i] = other.changed_mols[i];
        }
        
        ffcomponents = other.ffcomponents;
    }
    
    return *this;
}

/** Comparison operator - two forcefields are equal if they have the same
    UID and version number */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2GFF<Potential>::operator==(const Intra2B2GFF<Potential> &other) const
{
    return G2FF::operator==(other);
}

/** Comparison operator - two forcefields are equal if they have the same
    UID and version number */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2GFF<Potential>::operator!=(const Intra2B2GFF<Potential> &other) const
{
    return G2FF::operator!=(other);
}
    
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2B2GFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Intra2B2GFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2B2GFF<Potential>::what() const
{
    return Intra2B2GFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2B2GFF<Potential>* Intra2B2GFF<Potential>::clone() const
{
    return new Intra2B2GFF<Potential>(*this);
}

/** Assert that the intramolecular nonbonded scale factors for the 
    molecule with number 'molnum' are the same for the molecule
    in both groups in this forcefield
    
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::assertSameIntraScaleFactors(MolNum molnum) const
{
    if (mols[0].contains(molnum) and mols[1].contains(molnum))
    {
        const typename Potential::Molecule &mol0
                   = mols[0].moleculesByIndex()[mols[0].indexOf(molnum)];
                   
        const typename Potential::Molecule &mol1
                   = mols[1].moleculesByIndex()[mols[1].indexOf(molnum)];
        
        if (mol0.parameters().intraScaleFactors() != 
            mol1.parameters().intraScaleFactors())
        {
            SireFF::detail::throwIntra2B2GFFIncompatibeScaleFactorsError(
                    mol0.molecule().name(), mol0.molecule().number(), 
                    mol0.molecule().version(), mol1.molecule().version() ); 
        }
    }
}

/** Return the symbols representing the energy components of this forcefield */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const typename Potential::Components& Intra2B2GFF<Potential>::components() const
{
    return ffcomponents;
}

/** Function used to perform the work of changing the name of this 
    forcefield - this renames the component symbols and the molecule group */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::_pvt_updateName()
{
    ffcomponents = Components( this->name() );
    G2FF::_pvt_updateName();
}

/** Set the property called 'name' to the value 'value'

    \throw SireBase::missing_property
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2GFF<Potential>::setProperty(const QString &name, const Property &value)
{
    return Potential::setProperty(name, value);
}

/** Return the value of the property with name 'name'

    \throw SireBase::missing_property
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const Property& Intra2B2GFF<Potential>::property(const QString &name) const
{
    return Potential::property(name);
}

/** Return whether or not this forcefield contains a property
    called 'name' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2GFF<Potential>::containsProperty(const QString &name) const
{
    return Potential::containsProperty(name);
}

/** Return the properties available in this forcefield (and their values) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const Properties& Intra2B2GFF<Potential>::properties() const
{
    return Potential::properties();
}

/** Tell the forcefield that everything now should be calculated 
    from scratch (and thus to stop recording changes). This is a useful
    function to call if you are about to make many changes to this
    forcefield (e.g. performing an MD move), as this prevents the 
    potentially costly record keeping needed to calculate the change
    in energy associated with the move. */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::mustNowRecalculateFromScratch()
{
    //record that this forcefield is dirty
    G2FF::setDirty();
    
    //now clear any delta information
    for (int i=0; i<2; ++i)
    {
        changed_mols[i].clear();
    }
}

/** Return whether or not we need to record the changes to this   
    forcefield (not necessary if the energy has to be recalculated
    from scratch) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2GFF<Potential>::recordingChanges() const
{
    return not (G2FF::isDirty() and changed_mols[0].isEmpty() and
                                    changed_mols[1].isEmpty());
}

/** Record the change described in 'change' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::recordChange(quint32 groupid,
                    const typename Potential::ChangedMolecule &change)
{
    if (change.isEmpty())
        return;

    MolNum molnum = change.number();
    
    if (changed_mols[groupid].contains(molnum))
    {
        ChangedMolecule &old_change = changed_mols[groupid][molnum];
        
        if (old_change.oldMolecule() == change.newMolecule())
        {
            //we have reverted the change!
            changed_mols[groupid].remove(molnum);
            
            if (changed_mols[0].isEmpty() and changed_mols[1].isEmpty())    
                //there is now no change
                G2FF::setClean();
            
            return;
        }
        else
        {
            //this is yet another change
            changed_mols[groupid][molnum].change( change.newMolecule() );
        }
    }
    else
    {
        changed_mols[groupid].insert(molnum, change);
    }
    
    G2FF::setDirty();
}

/** Virtual function used to return the components of the forcefield
    (via the FF::components() base class function) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const FFComponent& Intra2B2GFF<Potential>::_pvt_components() const
{
    return this->components();
}

/** Virtual function called when the underlying potential energy surface
    has been changed (e.g. by changing the cutoff distance) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::changedPotential()
{
    G2FF::incrementVersion();
    this->mustNowRecalculateFromScratch();
}

/** Record the fact that the molecule 'mol' has been added to this forcefield 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::_pvt_added(quint32 groupid,
                                      const PartialMolecule &molecule, 
                                      const PropertyMap &map)
{
    Molecules old_state = mols[groupid];
    
    try
    {
        if (this->recordingChanges())
        {
            ChangedMolecule mol = mols[groupid].add(molecule, map, *this, true);
            this->assertSameIntraScaleFactors(molecule.number());
            this->recordChange(groupid, mol);
        }
        else
        {
            mols[groupid].add(molecule, map, *this, false);
            G2FF::setDirty();
            this->assertSameIntraScaleFactors(molecule.number());
        }
    }
    catch(...)
    {
        mols[groupid] = old_state;
        throw;
    }
}

/** Record the fact that the molecule 'mol' has been removed from this forcefield */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::_pvt_removed(quint32 groupid,
                                          const PartialMolecule &molecule)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols[groupid].remove(molecule, *this, true);
        this->recordChange(groupid, mol);
    }
    else
    {
        mols[groupid].remove(molecule, *this, false);
        G2FF::setDirty();
    }
}

/** Record that fact that the molecule 'molecule' has been changed

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::_pvt_changed(quint32 groupid,
                                          const SireMol::Molecule &molecule,
                                          bool auto_commit)
{
    Molecules old_state = mols[groupid];

    try
    {
        if (this->recordingChanges())
        {
            ChangedMolecule mol = mols[groupid].change(molecule, *this, false);
            this->assertSameIntraScaleFactors(molecule.number());
            
            this->recordChange(groupid, mol);
        }
        else
        {
            mols[groupid].change(molecule, *this, false);
            G2FF::setDirty();
            this->assertSameIntraScaleFactors(molecule.number());
        }
    }
    catch(...)
    {
        mols[groupid] = old_state;
        throw;
    }
}

/** Record that the provided list of molecules have changed 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::_pvt_changed(quint32 groupid,
                                          const QList<SireMol::Molecule> &molecules,
                                          bool auto_commit)
{
    Molecules old_mols = mols[groupid];
    QHash<MolNum,ChangedMolecule> old_changed_mols = changed_mols[groupid];

    try
    {
        if (this->recordingChanges())
        {   
            for (QList<SireMol::Molecule>::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                ChangedMolecule change = mols[groupid].change(*it, *this, true);
                this->assertSameIntraScaleFactors(it->number());
                
                this->recordChange(groupid, change);
            }
        }
        else
        {
            for (QList<SireMol::Molecule>::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                //we cannot add a molecule to this group that has a different
                //intramolecular non-bonded scale factor than that that already
                //exists...
                mols[groupid].change(*it, *this, false);
                this->assertSameIntraScaleFactors(it->number());
            }
            
            G2FF::setDirty();
        }
    }
    catch(...)
    {
        mols[groupid] = old_mols;
        changed_mols[groupid] = old_changed_mols;
        throw;
    }
}
    
/** Record that all of the molecules have been removed */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::_pvt_removedAll(quint32 groupid)
{
    mols[groupid].clear();

    this->mustNowRecalculateFromScratch();
}
 
/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */       
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2B2GFF<Potential>::_pvt_wouldChangeProperties(quint32 groupid,
                                                      MolNum molnum, 
                                                      const PropertyMap &map) const
{
    return mols[groupid].wouldChangeProperties(molnum, map);
}

/** Recalculate the energy of the current state of this forcefield. This
    will recalculate the energy using the quickest possible route, e.g.
    if will only recalculate the energies of molecules that have changed
    since the last evaluation */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2B2GFF<Potential>::recalculateEnergy()
{
    int nmols0 = mols[0].count();
    const ChunkedVector<typename Potential::Molecule> &mols0_array 
                            = mols[0].moleculesByIndex();

    int nmols1 = mols[1].count();
    const ChunkedVector<typename Potential::Molecule> &mols1_array
                            = mols[1].moleculesByIndex();

    if (changed_mols[0].count() == nmols0 or
        changed_mols[1].count() == nmols1)
    {
        //all of the molecules in a group have changed!
        changed_mols[0].clear();
        changed_mols[1].clear();
    }

    //tell the potential that we are starting an evaluation
    Potential::startEvaluation();

    try
    {

    if (changed_mols[0].isEmpty() and changed_mols[1].isEmpty())
    {
        //we are not recording changes, so we have to assume that
        //everything has changed. Recalculate the total energy from scratch
        EnergyWorkspace workspace;
        Energy total_nrg;

        //loop over all molecules in the smallest group...
        if (nmols0 <= nmols1)
        {
            for (int i=0; i<nmols0; ++i)
            {
                const typename Potential::Molecule &mol0 = mols0_array[i];
                
                //does this molecule exist in the the other group?
                if (not mols[1].contains(mol0.number()))
                    continue;
                    
                const typename Potential::Molecule &mol1
                                        = mols1_array[ mols[1].indexOf(mol0.number()) ];
                                        
                Potential::calculateEnergy(mol0, mol1, total_nrg, workspace);
            }
        }
        else
        {
            for (int j=0; j<nmols1; ++j)
            {
                const typename Potential::Molecule &mol1 = mols1_array[j];
                
                //does this molecule exist in the other group?
                if (not mols[0].contains(mol1.number()))
                    continue;
                    
                const typename Potential::Molecule &mol0
                                        = mols0_array[ mols[0].indexOf(mol1.number()) ];
                                        
                Potential::calculateEnergy(mol0, mol1, total_nrg, workspace);
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

        if (changed_mols[0].isEmpty())
        {
            for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                                it = changed_mols[1].constBegin();
                 it != changed_mols[1].constEnd();
                 ++it)
            {
                MolNum molnum = it.key();
                
                if (not mols[0].contains(molnum))
                    continue;
                    
                const typename Potential::Molecule &mol0 
                                        = mols0_array[ mols[0].indexOf(molnum) ];
                                        
                Potential::calculateEnergy(it->oldParts(), mol0,
                                           old_nrg, old_workspace);
                                           
                Potential::calculateEnergy(it->newParts(), mol0,
                                           new_nrg, new_workspace);
            }
        }
        else if (changed_mols[1].isEmpty())
        {
            for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                                it = changed_mols[0].constBegin();
                 it != changed_mols[0].constEnd();
                 ++it)
            {
                MolNum molnum = it.key();
                
                if (not mols[1].contains(molnum))
                    continue;
                    
                const typename Potential::Molecule &mol1 
                                        = mols1_array[ mols[1].indexOf(molnum) ];
                                        
                Potential::calculateEnergy(it->oldParts(), mol1,
                                           old_nrg, old_workspace);
                                           
                Potential::calculateEnergy(it->newParts(), mol1,
                                           new_nrg, new_workspace);
            }
        }
        else
        {
            //both changed!
            for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                                it = changed_mols[0].constBegin();
                 it != changed_mols[0].constEnd();
                 ++it)
            {
                MolNum molnum = it.key();
                
                typename QHash<MolNum,ChangedMolecule>::const_iterator
                                        it2 = changed_mols[1].constFind(molnum);
                
                if (it2 != changed_mols[1].constEnd())
                {
                    //both changed!
                    Potential::calculateEnergy(it->oldMolecule(),
                                               it2->oldMolecule(),
                                               old_nrg, old_workspace);
                                               
                    Potential::calculateEnergy(it->newMolecule(),
                                               it2->newMolecule(),
                                               new_nrg, new_workspace);
                }
                else if (mols[1].contains(molnum))
                {
                    const typename Potential::Molecule &mol1
                                         = mols1_array[ mols[1].indexOf(molnum) ];
                                         
                    Potential::calculateEnergy(it->oldParts(), mol1,
                                               old_nrg, old_workspace);
                    
                    Potential::calculateEnergy(it->newParts(), mol1,
                                               new_nrg, new_workspace);
                }
            }
            
            //there may be some molecules changed in mols[1] that we have missed...
            for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                   it = changed_mols[1].constBegin();
                 it != changed_mols[1].constEnd();
                 ++it)
            {
                MolNum molnum = it.key();
                
                if (not changed_mols[0].contains(molnum)
                    and mols[0].contains(molnum))
                {
                    const typename Potential::Molecule &mol0
                                     = mols0_array[ mols[0].indexOf(molnum) ];
                                     
                    Potential::calculateEnergy(mol0, it->oldParts(),
                                               old_nrg, old_workspace);
                                               
                    Potential::calculateEnergy(mol0, it->newParts(),
                                               new_nrg, new_workspace);
                }
            }
        }
         
        //change the energy
        this->components().changeEnergy(*this, new_nrg - old_nrg);
        
        //clear the changed molecules
        changed_mols[0].clear();
        changed_mols[1].clear();
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

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

namespace SireFF
{
namespace detail
{

template<class Potential>
struct Intra2B2GRMT
{
    static const RegisterMetaType< SireFF::Intra2B2GFF<Potential> > r_intra2b2gff;
};

template<class Potential>
const RegisterMetaType< SireFF::Intra2B2GFF<Potential> > 
Intra2B2GRMT<Potential>::r_intra2b2gff;

}
}

/** Serialise to a binary datastream */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const SireFF::Intra2B2GFF<Potential> &intra2b2gff)
{
    SireStream::writeHeader(ds, 
                SireFF::detail::Intra2B2GRMT<Potential>::r_intra2b2gff, 1);
    
    SireStream::SharedDataStream sds(ds);
    
    sds << intra2b2gff.mols[0] << intra2b2gff.mols[1] 
        << intra2b2gff.changed_mols[0] << intra2b2gff.changed_mols[1]
        << static_cast<const Potential&>(intra2b2gff)
        << static_cast<const SireFF::G2FF&>(intra2b2gff);
        
    return ds;
}

/** Extract from a binary datastream */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        SireFF::Intra2B2GFF<Potential> &intra2b2gff)
{
    SireStream::VersionID v = SireStream::readHeader(ds, 
                               SireFF::detail::Intra2B2GRMT<Potential>::r_intra2b2gff);
                                        
    if (v == 1)
    {
        SireStream::SharedDataStream sds(ds);
        
        sds >> intra2b2gff.mols[0] >> intra2b2gff.mols[1] 
            >> intra2b2gff.changed_mols[0] >> intra2b2gff.changed_mols[1]
            >> static_cast<Potential&>(intra2b2gff)
            >> static_cast<SireFF::G2FF&>(intra2b2gff);
            
        intra2b2gff._pvt_updateName();
        
        return ds;
    }
    else
        throw SireStream::version_error(v, "1",
                 SireFF::detail::Intra2B2GRMT<Potential>::r_intra2b2gff, CODELOC );

    return ds;
}

SIRE_END_HEADER

#endif
