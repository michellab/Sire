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

#ifndef SIREFF_INTER2BFF_HPP
#define SIREFF_INTER2BFF_HPP

#include "ff.h"
#include "g1ff.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireFF
{
template<class Potential>
class Inter2BFF;
}

template<class Potential>
QDataStream& operator<<(QDataStream&, const SireFF::Inter2BFF<Potential>&);

template<class Potential>
QDataStream& operator>>(QDataStream&, SireFF::Inter2BFF<Potential>&);

namespace SireFF
{

using SireBase::ChunkedVector;

/** This class provides an intermolecular non-bonded forcefield
    that can work with an two-body potential (provided via the 
    template type 'Potential'). 
    
    @author Christopher Woods
*/
template<class Potential>
class SIREFF_EXPORT Inter2BFF 
                : public SireBase::ConcreteProperty<Inter2BFF<Potential>, G1FF>, 
                  public Potential
{

friend QDataStream& ::operator<<<>(QDataStream&, const Inter2BFF<Potential>&);
friend QDataStream& ::operator>><>(QDataStream&, Inter2BFF<Potential>&);

public:
    typedef typename Potential::Components Components;

    Inter2BFF();
    Inter2BFF(const QString &name);
    
    Inter2BFF(const Inter2BFF<Potential> &other);
    
    ~Inter2BFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Inter2BFF<Potential>& operator=(const Inter2BFF<Potential> &other);
    
    bool operator==(const Inter2BFF<Potential> &other) const;
    bool operator!=(const Inter2BFF<Potential> &other) const;

    Inter2BFF<Potential>* clone() const;

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

    bool recordingChanges() const;

    void recordChange(const ChangedMolecule &change);

    void recalculateEnergy();

    void changedPotential();
    
    void _pvt_added(const PartialMolecule &mol, const PropertyMap &map);

    void _pvt_removed(const PartialMolecule &mol);

    void _pvt_changed(const SireMol::Molecule &molecule, bool auto_commit);
    void _pvt_changed(const QList<SireMol::Molecule> &molecules, bool auto_commit);
    
    void _pvt_removedAll();
        
    bool _pvt_wouldChangeProperties(MolNum molnum, 
                                    const PropertyMap &map) const;

    void _pvt_updateName();

    /** The parameterised version of the molecules in this forcefield */
    Molecules mols;

    /** The list of molecules that have changed since the last evaluation.
        While ffmols only contains the newest version of the molecule,
        this list contains both the newest version, and the version of the
        molecule at the last energy evaluation. */
    QHash<MolNum,ChangedMolecule> changed_mols;

    /** The numbers of the molecules that have been removed since
        the last energy evaluation */
    QSet<MolNum> removed_mols;

    /** The energy components available for this forcefield */
    Components ffcomponents;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2BFF<Potential>::Inter2BFF() 
                     : SireBase::ConcreteProperty<Inter2BFF<Potential>,G1FF>(), 
                       Potential()
{
    this->_pvt_updateName();
}

/** Construct this forcefield, providing it with a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2BFF<Potential>::Inter2BFF(const QString &ffname)
                     : SireBase::ConcreteProperty<Inter2BFF<Potential>,G1FF>(), 
                       Potential()
{
    FF::setName(ffname);
}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2BFF<Potential>::Inter2BFF(const Inter2BFF<Potential> &other)
                     : SireBase::ConcreteProperty<Inter2BFF<Potential>,G1FF>(other), 
                       Potential(other),
                       mols(other.mols), changed_mols(other.changed_mols),
                       removed_mols(other.removed_mols), 
                       ffcomponents(other.ffcomponents)
{}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2BFF<Potential>::~Inter2BFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2BFF<Potential>& 
Inter2BFF<Potential>::operator=(const Inter2BFF<Potential> &other)
{
    if (this != &other)
    {
        G1FF::operator=(other);
        Potential::operator=(other);
        
        mols = other.mols;
        changed_mols = other.changed_mols;
        removed_mols = other.removed_mols;
        ffcomponents = other.ffcomponents;
    }
    
    return *this;
}

/** Comparison operator - two forcefields are equal if they have the same
    UID and version number */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2BFF<Potential>::operator==(const Inter2BFF<Potential> &other) const
{
    return G1FF::operator==(other);
}

/** Comparison operator - two forcefields are equal if they have the same
    UID and version number */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2BFF<Potential>::operator!=(const Inter2BFF<Potential> &other) const
{
    return G1FF::operator!=(other);
}
    
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Inter2BFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Inter2BFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Inter2BFF<Potential>::what() const
{
    return Inter2BFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Inter2BFF<Potential>* Inter2BFF<Potential>::clone() const
{
    return new Inter2BFF<Potential>(*this);
}

/** Return the symbols representing the energy components of this forcefield */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const typename Potential::Components& Inter2BFF<Potential>::components() const
{
    return ffcomponents;
}

/** Function used to perform the work of changing the name of this 
    forcefield - this renames the component symbols and the molecule group */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::_pvt_updateName()
{
    ffcomponents = Components( this->name() );
    G1FF::_pvt_updateName();
}

/** Set the property called 'name' to the value 'value'

    \throw SireBase::missing_property
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2BFF<Potential>::setProperty(const QString &name, const Property &value)
{
    return Potential::setProperty(name, value);
}

/** Return the value of the property with name 'name'

    \throw SireBase::missing_property
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const Property& Inter2BFF<Potential>::property(const QString &name) const
{
    return Potential::property(name);
}

/** Return whether or not this forcefield contains a property
    called 'name' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2BFF<Potential>::containsProperty(const QString &name) const
{
    return Potential::containsProperty(name);
}

/** Return the properties available in this forcefield (and their values) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const Properties& Inter2BFF<Potential>::properties() const
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
void Inter2BFF<Potential>::mustNowRecalculateFromScratch()
{
    //record that this forcefield is dirty
    G1FF::setDirty();
    
    //now clear any delta information
    changed_mols.clear();
    removed_mols.clear();
}

/** Return whether or not we need to record the changes to this   
    forcefield (not necessary if the energy has to be recalculated
    from scratch) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2BFF<Potential>::recordingChanges() const
{
    return not (G1FF::isDirty() and changed_mols.isEmpty());
}

/** Record the change described in 'change' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::recordChange(
                    const typename Potential::ChangedMolecule &change)
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
            removed_mols.remove(molnum);
            
            if (changed_mols.isEmpty())
                //there are now no changes
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
    
    if (change.newMolecule().isEmpty())
        //the molecule has been removed
        removed_mols.insert(molnum);
    else
        //the molecule may have been re-added
        removed_mols.remove(molnum);

    G1FF::setDirty();
}

/** Virtual function called when the underlying potential energy surface
    has been changed (e.g. by changing the cutoff distance) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::changedPotential()
{
    G1FF::incrementVersion();
    this->mustNowRecalculateFromScratch();
}

/** Record the fact that the molecule 'mol' has been added to this forcefield 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::_pvt_added(const PartialMolecule &molecule, 
                                      const PropertyMap &map)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols.add(molecule, map, *this, true);
        this->recordChange(mol);
    }
    else
    {
        mols.add(molecule, map, *this, false);
        G1FF::setDirty();
    }
}

/** Record the fact that the molecule 'mol' has been removed from this forcefield */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::_pvt_removed(const PartialMolecule &molecule)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols.remove(molecule, *this, true);
        this->recordChange(mol);
    }
    else
    {
        mols.remove(molecule, *this, false);
        G1FF::setDirty();
    }
}

/** Record that fact that the molecule 'molecule' has been changed

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::_pvt_changed(const SireMol::Molecule &molecule, bool auto_commit)
{
    if (this->recordingChanges())
    {
        ChangedMolecule mol = mols.change(molecule, *this, true);
        this->recordChange(mol);
    }
    else
    {
        mols.change(molecule, *this, false);
        G1FF::setDirty();
    }
}

/** Record that the provided list of molecules have changed 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::_pvt_changed(const QList<SireMol::Molecule> &molecules,
                                        bool auto_commit)
{
    Molecules old_mols = mols;
    QHash<MolNum,ChangedMolecule> old_changed_mols = changed_mols;
    QSet<MolNum> old_removed_mols = removed_mols;

    try
    {
        if (this->recordingChanges())
        {   
            for (QList<SireMol::Molecule>::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                ChangedMolecule change = mols.change(*it, *this, true);
                this->recordChange(change);
            }
        }
        else
        {
            for (QList<SireMol::Molecule>::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                mols.change(*it, *this, false);
            }

            G1FF::setDirty();
        }
    }
    catch(...)
    {
        mols = old_mols;
        changed_mols = old_changed_mols;
        removed_mols = old_removed_mols;
        throw;
    }
}
    
/** Record that all of the molecules have been removed */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::_pvt_removedAll()
{
    mols.clear();
    this->mustNowRecalculateFromScratch();
}
 
/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */       
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Inter2BFF<Potential>::_pvt_wouldChangeProperties(MolNum molnum, 
                                                      const PropertyMap &map) const
{
    return mols.wouldChangeProperties(molnum, map);
}

/** Recalculate the energy of the current state of this forcefield. This
    will recalculate the energy using the quickest possible route, e.g.
    if will only recalculate the energies of molecules that have changed
    since the last evaluation */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Inter2BFF<Potential>::recalculateEnergy()
{
    int nmols = mols.count();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                                                = mols.moleculesByIndex();

    if (changed_mols.count() == nmols)
        //all of the molecules have changed!
        changed_mols.clear();

    //tell the potential that we are starting an evaluation
    Potential::startEvaluation();

    try
    {

    if (changed_mols.isEmpty())
    {
        Energy total_nrg;

        {
            //we are not recording changes, so we have to assume that
            //everything has changed. Recalculate the total energy from scratch
            EnergyWorkspace workspace;
            Energy my_total_nrg;

            const ChunkedVector<typename Potential::Molecule> &my_mols_array = mols_array;
            const int my_nmols = nmols;

            //loop over all pairs of molecules
            for (int i=0; i<my_nmols-1; ++i)
            {
                const typename Potential::Molecule &mol0 = my_mols_array[i];
        
                for (int j=i+1; j<my_nmols; ++j)
                {
                    const typename Potential::Molecule &mol1 = my_mols_array[j];
                    Potential::calculateEnergy(mol0, mol1, my_total_nrg, workspace);
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

        EnergyWorkspace old_workspace;
        EnergyWorkspace new_workspace;
        
        for (int i=0; i<nmols; ++i)
        {
            const typename Potential::Molecule &mol0 = mols_array[i];
            
            typename QHash<MolNum,ChangedMolecule>::const_iterator it
                                               = changed_mols.constFind(mol0.number());
                                               
            if (it == changed_mols.constEnd())
            {
                //this molecule has not changed - just calculate its
                //energy with all of the changed molecules - as this molecule
                //hasn't changed, we only need to calculate the change
                //in energy between this molecule and the changed parts
                //of the changed molecules
                for (typename QHash<MolNum,ChangedMolecule>::const_iterator
                                      it2 = changed_mols.constBegin();
                     it2 != changed_mols.constEnd();
                     ++it2)
                {
                    Potential::calculateEnergy(mol0, it2->oldParts(),
                                               old_nrg, old_workspace);
                                               
                    Potential::calculateEnergy(mol0, it2->newParts(),
                                               new_nrg, new_workspace);
                }
            }
            else if (changed_mols.count() > 1)
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
                    for (++it2; it2 != changed_mols.constEnd(); ++it2)
                    {
                        Potential::calculateEnergy(it->oldMolecule(),
                                                   it2->oldMolecule(),
                                                   old_nrg, old_workspace);
                                                   
                        Potential::calculateEnergy(it->newMolecule(),
                                                   it2->newMolecule(),
                                                   new_nrg, new_workspace);
                    }
                }
                else
                {
                    for (++it2; it2 != changed_mols.constEnd(); ++it2)
                    {
                        if (it2->changedAll())
                        {
                            //all of the other molecule has changed - we need
                            //to calculate the energy of the whole molecules
                            //interaction
                            Potential::calculateEnergy(it->oldMolecule(),
                                                       it2->oldMolecule(),
                                                       old_nrg, old_workspace);
                            
                            Potential::calculateEnergy(it->newMolecule(),
                                                       it2->newMolecule(),
                                                       new_nrg, new_workspace);
                        }
                        else
                        {
                            //both a part of this molecule and a part of the 
                            //other molecule have changed
                           
                            //the change in energy associated with changing 
                            //the first molecule...
                            Potential::calculateEnergy(it->oldParts(),
                                                       it2->oldMolecule(),
                                                       old_nrg, old_workspace);
                                                       
                            Potential::calculateEnergy(it->newParts(),
                                                       it2->newMolecule(),
                                                       new_nrg, new_workspace);
                                                       
                           //now the change in energy associated with changing
                           //the second molecule...
                           Potential::calculateEnergy(it2->oldParts(),
                                                      it->oldMolecule(),
                                                      old_nrg, old_workspace);
                                                     
                           Potential::calculateEnergy(it2->newParts(),
                                                      it->newMolecule(),
                                                      new_nrg, new_workspace);
                                                      
                           //now remove double counted changed in mol1 with
                           //change in mol2
                           Potential::calculateEnergy(it->oldParts(),
                                                      it2->oldParts(),
                                                      old_nrg, old_workspace, -1);
                                                      
                           Potential::calculateEnergy(it->newParts(),
                                                      it2->newParts(),
                                                      new_nrg, new_workspace, -1);
                        }
                    }
                }
            }
        }

        //finally, loop over all of the molecules that have been removed - the energy
        //of non-changed molecules with removed molecules has already been calculated,
        //as has the energy of moved molecules that are before the removed molecules
        //in the moved list. We only now have to calculate the energy of the removed
        //molecules with all of the molecules that lie above us in the moved list
        for (typename QSet<MolNum>::const_iterator it = removed_mols.constBegin();
             it != removed_mols.constEnd();
             ++it)
        {
            typename QSet<MolNum>::const_iterator it2 = it;
            
            const ChangedMolecule &mol0 = *(changed_mols.constFind(*it));
            
            for (++it2; it2 != removed_mols.constEnd(); ++it2)
            {
                const ChangedMolecule &mol1 = *(changed_mols.constFind(*it2));
            
                Potential::calculateEnergy(mol0.oldMolecule(),
                                           mol1.oldMolecule(),
                                           old_nrg, old_workspace);
                                           
                //molecule has been removed, so no new energy
            }
        }
         
        //change the energy
        this->components().changeEnergy(*this, new_nrg - old_nrg);
        
        //clear the changed molecules
        changed_mols.clear();
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
struct Inter2BRMT
{
    static const RegisterMetaType< SireFF::Inter2BFF<Potential> > r_inter2bff;
};

template<class Potential>
const RegisterMetaType< SireFF::Inter2BFF<Potential> > 
Inter2BRMT<Potential>::r_inter2bff;

}
}

/** Serialise to a binary datastream */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const SireFF::Inter2BFF<Potential> &inter2bff)
{
    SireStream::writeHeader(ds, SireFF::detail::Inter2BRMT<Potential>::r_inter2bff, 1);
    
    SireStream::SharedDataStream sds(ds);
    
    sds << inter2bff.mols << inter2bff.changed_mols << inter2bff.removed_mols
        << static_cast<const Potential&>(inter2bff)
        << static_cast<const SireFF::G1FF&>(inter2bff);
        
    return ds;
}

/** Extract from a binary datastream */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        SireFF::Inter2BFF<Potential> &inter2bff)
{
    SireStream::VersionID v = SireStream::readHeader(ds, 
                                 SireFF::detail::Inter2BRMT<Potential>::r_inter2bff);
                                        
    if (v == 1)
    {
        SireStream::SharedDataStream sds(ds);
        
        sds >> inter2bff.mols >> inter2bff.changed_mols >> inter2bff.removed_mols
            >> static_cast<Potential&>(inter2bff)
            >> static_cast<SireFF::G1FF&>(inter2bff);
            
        inter2bff._pvt_updateName();
        
        return ds;
    }
    else
        throw SireStream::version_error(v, "1",
                     SireFF::detail::Inter2BRMT<Potential>::r_inter2bff, CODELOC );

    return ds;
}

SIRE_END_HEADER

#endif
