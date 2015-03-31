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

#ifndef SIREFF_INTRA2BFF_HPP
#define SIREFF_INTRA2BFF_HPP

#include "ff.h"
#include "g1ff.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireFF
{
template<class Potential>
class Intra2BFF;
}

template<class Potential>
QDataStream& operator<<(QDataStream&, const SireFF::Intra2BFF<Potential>&);

template<class Potential>
QDataStream& operator>>(QDataStream&, SireFF::Intra2BFF<Potential>&);

namespace SireFF
{

using SireBase::ChunkedVector;

/** This class provides an intramolecular non-bonded forcefield
    that can work with an two-body potential (provided via the 
    template type 'Potential'). 
    
    @author Christopher Woods
*/
template<class Potential>
class SIREFF_EXPORT Intra2BFF 
                : public SireBase::ConcreteProperty<Intra2BFF<Potential>, G1FF>, 
                  public Potential
{

friend QDataStream& ::operator<<<>(QDataStream&, const Intra2BFF<Potential>&);
friend QDataStream& ::operator>><>(QDataStream&, Intra2BFF<Potential>&);

public:
    typedef typename Potential::Components Components;

    Intra2BFF();
    Intra2BFF(const QString &name);
    
    Intra2BFF(const Intra2BFF<Potential> &other);
    
    ~Intra2BFF();
    
    static const char* typeName();
    
    const char* what() const;
    
    Intra2BFF<Potential>& operator=(const Intra2BFF<Potential> &other);
    
    bool operator==(const Intra2BFF<Potential> &other) const;
    bool operator!=(const Intra2BFF<Potential> &other) const;

    Intra2BFF<Potential>* clone() const;

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

    const FFComponent& _pvt_components() const;

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

    /** The energy components available for this forcefield */
    Components ffcomponents;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2BFF<Potential>::Intra2BFF() 
                     : SireBase::ConcreteProperty<Intra2BFF<Potential>,G1FF>(), 
                       Potential()
{
    this->_pvt_updateName();
}

/** Construct this forcefield, providing it with a name */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2BFF<Potential>::Intra2BFF(const QString &ffname)
                     : SireBase::ConcreteProperty<Intra2BFF<Potential>,G1FF>(), 
                       Potential()
{
    FF::setName(ffname);
}

/** Copy constructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2BFF<Potential>::Intra2BFF(const Intra2BFF<Potential> &other)
                     : SireBase::ConcreteProperty<Intra2BFF<Potential>,G1FF>(other), 
                       Potential(other),
                       mols(other.mols), changed_mols(other.changed_mols),
                       ffcomponents(other.ffcomponents)
{}

/** Destructor */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2BFF<Potential>::~Intra2BFF()
{}

/** Copy assignment operator */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2BFF<Potential>& 
Intra2BFF<Potential>::operator=(const Intra2BFF<Potential> &other)
{
    if (this != &other)
    {
        G1FF::operator=(other);
        Potential::operator=(other);
        
        mols = other.mols;
        changed_mols = other.changed_mols;
        ffcomponents = other.ffcomponents;
    }
    
    return *this;
}

/** Comparison operator - two forcefields are equal if they have the same
    UID and version number */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2BFF<Potential>::operator==(const Intra2BFF<Potential> &other) const
{
    return G1FF::operator==(other);
}

/** Comparison operator - two forcefields are equal if they have the same
    UID and version number */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2BFF<Potential>::operator!=(const Intra2BFF<Potential> &other) const
{
    return G1FF::operator!=(other);
}
    
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2BFF<Potential>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Intra2BFF<Potential> >() );
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const char* Intra2BFF<Potential>::what() const
{
    return Intra2BFF<Potential>::typeName();
}

template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
Intra2BFF<Potential>* Intra2BFF<Potential>::clone() const
{
    return new Intra2BFF<Potential>(*this);
}

/** Return the symbols representing the energy components of this forcefield */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const typename Potential::Components& Intra2BFF<Potential>::components() const
{
    return ffcomponents;
}

/** Function used to perform the work of changing the name of this 
    forcefield - this renames the component symbols and the molecule group */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2BFF<Potential>::_pvt_updateName()
{
    ffcomponents = Components( this->name() );
    G1FF::_pvt_updateName();
}

/** Set the property called 'name' to the value 'value'

    \throw SireBase::missing_property
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2BFF<Potential>::setProperty(const QString &name, const Property &value)
{
    return Potential::setProperty(name, value);
}

/** Return the value of the property with name 'name'

    \throw SireBase::missing_property
*/
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const Property& Intra2BFF<Potential>::property(const QString &name) const
{
    return Potential::property(name);
}

/** Return whether or not this forcefield contains a property
    called 'name' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2BFF<Potential>::containsProperty(const QString &name) const
{
    return Potential::containsProperty(name);
}

/** Return the properties available in this forcefield (and their values) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const Properties& Intra2BFF<Potential>::properties() const
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
void Intra2BFF<Potential>::mustNowRecalculateFromScratch()
{
    //record that this forcefield is dirty
    G1FF::setDirty();
    
    //now clear any delta information
    changed_mols.clear();
}

/** Return whether or not we need to record the changes to this   
    forcefield (not necessary if the energy has to be recalculated
    from scratch) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2BFF<Potential>::recordingChanges() const
{
    return not (G1FF::isDirty() and changed_mols.isEmpty());
}

/** Record the change described in 'change' */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2BFF<Potential>::recordChange(
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
            
            if (changed_mols.isEmpty())
                //there are no changes
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

/** Virtual function used to return the components of the forcefield
    (via the FF::components() base class function) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
const FFComponent& Intra2BFF<Potential>::_pvt_components() const
{
    return this->components();
}

/** Virtual function called when the underlying potential energy surface
    has been changed (e.g. by changing the cutoff distance) */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2BFF<Potential>::changedPotential()
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
void Intra2BFF<Potential>::_pvt_added(const PartialMolecule &molecule, 
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
void Intra2BFF<Potential>::_pvt_removed(const PartialMolecule &molecule)
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
void Intra2BFF<Potential>::_pvt_changed(const SireMol::Molecule &molecule, bool auto_commit)
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
void Intra2BFF<Potential>::_pvt_changed(const QList<SireMol::Molecule> &molecules,
                                        bool auto_commit)
{
    Molecules old_mols = mols;
    QHash<MolNum,ChangedMolecule> old_changed_mols = changed_mols;

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
        throw;
    }
}
    
/** Record that all of the molecules have been removed */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
void Intra2BFF<Potential>::_pvt_removedAll()
{
    mols.clear();
    this->mustNowRecalculateFromScratch();
}
 
/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */       
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
bool Intra2BFF<Potential>::_pvt_wouldChangeProperties(MolNum molnum, 
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
void Intra2BFF<Potential>::recalculateEnergy()
{
    int nmols = mols.count();
    const ChunkedVector<typename Potential::Molecule> &mols_array 
                            = mols.moleculesByIndex();

    //tell the potential that we are starting an evaluation
    Potential::startEvaluation();

    try
    {

    if (changed_mols.isEmpty())
    {
        //we are not recording changes, so we have to assume that
        //everything has changed. Recalculate the total energy from scratch
        EnergyWorkspace workspace;
        Energy total_nrg;

        //loop over all molecules and calculate their intramolecular energies
        for (int i=0; i<nmols; ++i)
        {
            const typename Potential::Molecule &mol = mols_array[i];
            Potential::calculateEnergy(mol, total_nrg, workspace);
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

        for (typename QHash<MolNum,ChangedMolecule>::const_iterator it
                                                = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (it->changedAll())
            {
                //calculate the change in energy of this whole molecule
                Potential::calculateEnergy(it->oldMolecule(),
                                           old_nrg, old_workspace);
                        
                Potential::calculateEnergy(it->newMolecule(),
                                           new_nrg, new_workspace);
            }
            else
            {
                //calculate the change in energy of the changed part
                //of the molecule
                Potential::calculateEnergy(it->oldParts(), old_nrg, old_workspace);

                Potential::calculateEnergy(it->oldParts(), it->oldMolecule(),
                                           old_nrg, old_workspace);
                                           
                Potential::calculateEnergy(it->newParts(), new_nrg, new_workspace);
                                           
                Potential::calculateEnergy(it->newParts(), it->newMolecule(),
                                           new_nrg, new_workspace);
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
struct Intra2BRMT
{
    static const RegisterMetaType< SireFF::Intra2BFF<Potential> > r_intra2bff;
};

template<class Potential>
const RegisterMetaType< SireFF::Intra2BFF<Potential> > 
Intra2BRMT<Potential>::r_intra2bff;

}
}

/** Serialise to a binary datastream */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const SireFF::Intra2BFF<Potential> &intra2bff)
{
    SireStream::writeHeader(ds, SireFF::detail::Intra2BRMT<Potential>::r_intra2bff, 1);
    
    SireStream::SharedDataStream sds(ds);
    
    sds << intra2bff.mols << intra2bff.changed_mols
        << static_cast<const Potential&>(intra2bff)
        << static_cast<const SireFF::G1FF&>(intra2bff);
        
    return ds;
}

/** Extract from a binary datastream */
template<class Potential>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        SireFF::Intra2BFF<Potential> &intra2bff)
{
    SireStream::VersionID v = SireStream::readHeader(ds, 
                                SireFF::detail::Intra2BRMT<Potential>::r_intra2bff);
                                        
    if (v == 1)
    {
        SireStream::SharedDataStream sds(ds);
        
        sds >> intra2bff.mols >> intra2bff.changed_mols
            >> static_cast<Potential&>(intra2bff)
            >> static_cast<SireFF::G1FF&>(intra2bff);
            
        intra2bff._pvt_updateName();
        
        return ds;
    }
    else
        throw SireStream::version_error(v, "1",
                    SireFF::detail::Intra2BRMT<Potential>::r_intra2bff, CODELOC );

    return ds;
}

SIRE_END_HEADER

#endif
