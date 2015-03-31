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

#ifndef SIREFF_DETAIL_FFMOLECULES_H
#define SIREFF_DETAIL_FFMOLECULES_H

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"

#include "SireBase/propertymap.h"
#include "SireBase/chunkedvector.hpp"

#include "SireStream/shareddatastream.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireFF
{

namespace detail
{
class FFMoleculeBase;
class FFMoleculesBase;

template<class PTNL>
class FFMolecule;

template<class PTNL>
class FFMolecules;

template<class FFMOL>
class ChangedMolecule;
}

}

QDataStream& operator<<(QDataStream&, const SireFF::detail::FFMoleculeBase&);
QDataStream& operator>>(QDataStream&, SireFF::detail::FFMoleculeBase&);

QDataStream& operator<<(QDataStream&, const SireFF::detail::FFMoleculesBase&);
QDataStream& operator>>(QDataStream&, SireFF::detail::FFMoleculesBase&);

template<class PTNL>
QDataStream& operator<<(QDataStream&, const SireFF::detail::FFMolecule<PTNL>&);
template<class PTNL>
QDataStream& operator>>(QDataStream&, SireFF::detail::FFMolecule<PTNL>&);

template<class PTNL>
QDataStream& operator<<(QDataStream&, const SireFF::detail::FFMolecules<PTNL>&);
template<class PTNL>
QDataStream& operator>>(QDataStream&, SireFF::detail::FFMolecules<PTNL>&);

template<class FFMOL>
QDataStream& operator<<(QDataStream&, const SireFF::detail::ChangedMolecule<FFMOL>&);
template<class FFMOL>
QDataStream& operator>>(QDataStream&, SireFF::detail::ChangedMolecule<FFMOL>&);

namespace SireMol
{
class MolNum;
}

namespace SireFF
{

namespace detail
{

using SireMol::MoleculeGroup;
using SireMol::PartialMolecule;
using SireMol::AtomSelection;
using SireMol::MolNum;
using SireMol::CGIdx;

using SireBase::PropertyMap;
using SireBase::ChunkedVector;

void throwFFMoleculesIncompatibleParameterNames(const PropertyMap &old_map, 
                                                const PropertyMap &new_map);

void assertCompatible(const FFMoleculeBase &mol0, const FFMoleculeBase &mol1);

/** This is the base class of most (all?) of the forcefield
    specialised molecules. Forcefields used specialised molecule
    classes as a way of storing important data in a way that
    is best for the forcefield itself. As these specialised
    classes are only used within a forcefield, they do not
    form part of the public API of Sire, and should not be
    used in any other code.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFMoleculeBase
{

friend QDataStream& ::operator<<(QDataStream&, const FFMoleculeBase&);
friend QDataStream& ::operator>>(QDataStream&, FFMoleculeBase&);

public:
    FFMoleculeBase();
    FFMoleculeBase(const PartialMolecule &molview);
    
    ~FFMoleculeBase();
    
    bool operator==(const FFMoleculeBase &other) const;
    bool operator!=(const FFMoleculeBase &other) const;
    
    int nCutGroups() const;
    
    CGIdx cgIdx(quint32 i) const;
    
    const PartialMolecule& molecule() const;

    MolNum number() const;

    bool isEmpty() const;

protected:
    FFMoleculeBase& operator=(const FFMoleculeBase &other);

    bool change(const SireMol::Molecule &molecule);
    bool add(const AtomSelection &added_atoms);
    bool remove(const AtomSelection &removed_atoms);
    
    void restore(const PartialMolecule &oldmol);
    
private:
    /** Copy of the view that this object represents */
    PartialMolecule mol;
    
    /** List of the CutGroups that are selected as part of this 
        molecule, in the order that they appear */
    QList<CGIdx> idx_to_cgidx;
};

/** This is the forcefield specific part of the FFMolecule class
    that holds the forcefield parameters for this part of the molecule
    
    @author Christopher Woods
*/
template<class PTNL>
class FFMolecule : public FFMoleculeBase
{

friend QDataStream& ::operator<<<>(QDataStream&, const FFMolecule<PTNL>&);
friend QDataStream& ::operator>><>(QDataStream&, FFMolecule<PTNL>&);

public:
    typedef typename PTNL::Parameters Parameters;
    typedef typename PTNL::ParameterNames ParameterNames;

    FFMolecule();
    
    FFMolecule(const PartialMolecule &molecule,
               PTNL &forcefield,
               const PropertyMap &map = PropertyMap());

    FFMolecule(const FFMolecule<PTNL> &other);

    ~FFMolecule();

    FFMolecule<PTNL>& operator=(const FFMolecule<PTNL> &other);
    
    bool operator==(const FFMolecule<PTNL> &other) const;
    bool operator!=(const FFMolecule<PTNL> &other) const;

    const Parameters& parameters() const;
    
    bool change(const SireMol::Molecule &molecule,
                PTNL &forcefield,
                const PropertyMap &map);

    bool add(const AtomSelection &selected_atoms,
             PTNL &forcefield,
             const PropertyMap &map);

    bool remove(const AtomSelection &selected_atoms,
                PTNL &forcefield,
                const PropertyMap &map);

    FFMolecule<PTNL> getDifferences(const FFMolecule<PTNL> &other) const;

protected:
    FFMolecule(const PartialMolecule &molecule,
               const Parameters &parameters);

    Parameters& _edit_parameters();

private:
    /** The parameters for this molecule */
    Parameters params;
};

/** This is the base class of most (all?) of the 'Molecules'
    classes. These classes hold groups of molecules together
    in a way that allows very rapid iteration over them and 
    their forcefield parameters. By its nature, this is a 
    template class, so this provides the template independent
    parts.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFMoleculesBase
{

friend QDataStream& ::operator<<(QDataStream&, const FFMoleculesBase&);
friend QDataStream& ::operator>>(QDataStream&, FFMoleculesBase&);

public:
    ~FFMoleculesBase();

    int count() const;

    bool isEmpty() const;

    const QVector<MolNum> molNumsByIndex() const;
    const QHash<MolNum,quint32>& indexesByMolNum() const;

    const QVector<PropertyMap>& parameterNamesByIndex() const;

    quint32 indexOf(MolNum molnum) const;
    bool contains(MolNum molnum) const;

    void assertContains(MolNum molnum) const;

protected:
    FFMoleculesBase();
    FFMoleculesBase(const MoleculeGroup &molgroup, const PropertyMap &map);
    
    FFMoleculesBase(const FFMoleculesBase &other);
    
    FFMoleculesBase& operator=(const FFMoleculesBase &other);

    bool operator==(const FFMoleculesBase &other) const;
    bool operator!=(const FFMoleculesBase &other) const;

    quint32 _pvt_add(MolNum molnum);

    void _pvt_remove(MolNum molnum);
    void _pvt_remove(const QList<MolNum> &molnums);

    void clear();

    /** The names of the properties used to get the parameters for
        each molecule */
    QVector<PropertyMap> parameter_names;
    
private:
    void _pvt_reindex();

    /** The number of each molecule in the group, 
        in the order they appear in the arrays */
    QVector<MolNum> molnums_by_idx;
    
    /** The indicies of each molecule in the array
        indexed by molecule number */
    QHash<MolNum,quint32> idxs_by_molnum;
};

/** This implementation class holds the information about 
    a change in a molecule. This can be used by a
    forcefield to simplify the calculation of the change in
    energy associated with a change in the molecule(s)
    (e.g. if only changes in energy associated with the 
    changed molecules needs to be evaluated) 
    
    The template parameter 'FFMOL' is the specialised Molecule
    type used by the forcefield. 
    
    @author Christopher Woods
*/
template<class FFMOL>
class ChangedMolecule
{

friend QDataStream& ::operator<<<>(QDataStream&, const ChangedMolecule<FFMOL>&);
friend QDataStream& ::operator>><>(QDataStream&, ChangedMolecule<FFMOL>&);

public:
    typedef FFMOL Molecule;

    ChangedMolecule();
    
    ChangedMolecule(const FFMOL &molecule);
    ChangedMolecule(const FFMOL &oldmol, const FFMOL &newmol);
    
    ChangedMolecule(const ChangedMolecule<FFMOL> &other);
    
    ~ChangedMolecule();
    
    ChangedMolecule<FFMOL>& operator=(const ChangedMolecule<FFMOL> &other);
    
    bool operator==(const ChangedMolecule<FFMOL> &other) const;
    bool operator!=(const ChangedMolecule<FFMOL> &other) const;
    
    MolNum number() const;
    
    bool isEmpty() const;
    
    bool changedAll() const;
    bool nothingChanged() const;
    
    const FFMOL& oldMolecule() const;
    const FFMOL& newMolecule() const;
    
    const FFMOL& oldParts() const;
    const FFMOL& newParts() const;
    
    ChangedMolecule<FFMOL>& removed();
    
    ChangedMolecule<FFMOL>& change(const FFMOL &new_molecule);
    
    static ChangedMolecule<FFMOL> removeMolecule(const FFMOL &molecule);
    
private:
    /** The complete, old molecule */
    FFMOL old_molecule;
    
    /** The complete, new molecule */
    FFMOL new_molecule;
    
    /** The parts of the old molecule that have changed */
    FFMOL old_parts;
    
    /** The corresponding changed parts of the new molecule */
    FFMOL new_parts;
};

/** This class provides an array of FF specialised molecules, 
    arranged in such a way as to speed up indexing over them.
    
    The template parameter is the potential type for which 
    this group is used (e.g. SireMM::CLJPotential)
    
    This is necessary, as different
    potentials use different parameterisation methods,
    and can supply different specialised molecule types.
    
    @author Christopher Woods
*/
template<class PTNL>
class FFMolecules : public FFMoleculesBase
{

friend QDataStream& ::operator<<<>(QDataStream&, const FFMolecules<PTNL>&);
friend QDataStream& ::operator>><>(QDataStream&, FFMolecules<PTNL>&);

public:
    typedef typename PTNL::Molecule Molecule;
    typedef SireFF::detail::ChangedMolecule<Molecule> ChangedMolecule;

    typedef typename Molecule::Parameters Parameters;
    typedef typename Molecule::ParameterNames ParameterNames;

    FFMolecules();
    
    FFMolecules(const MoleculeGroup &molgroup, PTNL &forcefield,
                const PropertyMap &map = PropertyMap());
    
    FFMolecules(const FFMolecules<PTNL> &other);
    
    ~FFMolecules();

    FFMolecules<PTNL>& operator=(const FFMolecules<PTNL> &other);
    
    bool operator==(const FFMolecules<PTNL> &other) const;
    bool operator!=(const FFMolecules<PTNL> &other) const;

    ChangedMolecule change(const SireMol::Molecule &molecule,
                           PTNL &forcefield,
                           bool record_changes = true);
                           
    ChangedMolecule add(const PartialMolecule &molecule,
                        const PropertyMap &map,
                        PTNL &forcefield,
                        bool record_changes = true);
                        
    ChangedMolecule remove(const PartialMolecule &molecule,
                           PTNL &forcefield,
                           bool record_changes = true);
    
    bool wouldChangeProperties(MolNum molnum, const PropertyMap &map) const;
    
    const ChunkedVector<Molecule>& moleculesByIndex() const;
    
    void clear();
    
protected:
    /** The array of forcefield-specialised molecules in this group */
    ChunkedVector<Molecule> mols_by_idx;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

////////
//////// Implementation of FFMolecule
////////

/** Null constructor */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL>::FFMolecule() : FFMoleculeBase()
{}

/** Create a specialised version of the molecule 'molecule'
    that will be part of the passed forcefield 'forcefield',
    using the parameter names in 'parameternames' to find
    and extract the correct parameters */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL>::FFMolecule(const PartialMolecule &molecule,
                             PTNL &forcefield,
                             const PropertyMap &map)
                 : FFMoleculeBase(molecule),
                   params( forcefield.getParameters(molecule,map) )
{}

/** Private constructor used to create a molecule from the passed
    data and parameters - these are assumed to be compatible - this
    is not checked! */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL>::FFMolecule(const PartialMolecule &molecule,
                             const typename FFMolecule<PTNL>::Parameters &parameters)
                 : FFMoleculeBase(molecule),
                   params(parameters)
{}

/** Copy constructor */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL>::FFMolecule(const FFMolecule<PTNL> &other)
                 : FFMoleculeBase(other), params(other.params)
{}

/** Destructor */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL>::~FFMolecule()
{}

/** Copy assignment operator */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL>& FFMolecule<PTNL>::operator=(const FFMolecule<PTNL> &other)
{
    FFMoleculeBase::operator=(other);
    params = other.params;
    
    return *this;
}

/** Comparison operator */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecule<PTNL>::operator==(const FFMolecule<PTNL> &other) const
{
    return FFMoleculeBase::operator==(other) and 
           params == other.params;
}

/** Comparison operator */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecule<PTNL>::operator!=(const FFMolecule<PTNL> &other) const
{
    return FFMoleculeBase::operator!=(other) or
           params != other.params;
}

/** Return the parameters for this view of the molecule. These parameters
    are specialised for the particular forcefield that this object
    was parameterised in (so will not necessarily make sense outside
    that forcefield) */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
const typename FFMolecule<PTNL>::Parameters& FFMolecule<PTNL>::parameters() const
{
    return params;
}

/** Protected function used by derived classes to edit the parameters */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
typename FFMolecule<PTNL>::Parameters& FFMolecule<PTNL>::_edit_parameters()
{
    return params;
}

/** Change this molecule so that it is equal to 'new_molecule'.

    Note that you can only change a molecule's molecule layout ID
    if *all* of the molecule is contained in this forcefield.
    If only part of the molecule is contained, then you will
    need to remove it, and then add it again to be able to
    change its layout ID.

    The names of the properties used to originally get the 
    parameters for this molecule are in 'parameternames'
    
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
    \throw SireBase::missing_property
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecule<PTNL>::change(const SireMol::Molecule &new_molecule,
                              PTNL &forcefield, const PropertyMap &old_map)
{
    //save the old molecule in case it needs to be restored
    PartialMolecule old_molecule = FFMoleculeBase::molecule();
    
    //change the molecule itself...
    bool mol_changed = FFMoleculeBase::change(new_molecule);
    
    if (mol_changed)
    {
        try
        {
            //now try to change the parameters
            params = forcefield.updateParameters(params, old_molecule,
                                                 FFMoleculeBase::molecule(),
                                                 old_map);
        }
        catch(...)
        {
            //something went wrong... restore the old molecule
            FFMoleculeBase::restore(old_molecule);
            throw;
        }
        
        return true;
    }
    else
        return false;
}

/** Add the selected atoms in 'added_atoms' to this view.

    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
    \throw SireBase::missing_property
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecule<PTNL>::add(const AtomSelection &added_atoms,
                           PTNL &forcefield, const PropertyMap &map)
{
    //save the old version of the molecule
    PartialMolecule old_mol = FFMoleculeBase::molecule();
    
    //add the atoms
    if (FFMoleculeBase::add(added_atoms))
    {
        //atoms were added. We now need to reparameterise the molecule
        try
        {
            params = forcefield.getParameters(FFMoleculeBase::molecule(),
                                              map);
        }
        catch(...)
        {
            FFMoleculeBase::restore(old_mol);
            throw;
        }
        
        return true;
    }
    else
        return false;
}

/** Remove the selected atoms 'removed_atoms' from this molecule.

    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
    \throw SireBase::missing_property
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecule<PTNL>::remove(const AtomSelection &removed_atoms,
                              PTNL &forcefield, const PropertyMap &map)
{
    if (this->isEmpty())
        return false;

    PartialMolecule old_mol = FFMoleculeBase::molecule();
    
    if (FFMoleculeBase::remove(removed_atoms))
    {
        if (FFMoleculeBase::molecule().selection().selectedNone())
        {
            //the molecule has been completely removed
            params = Parameters();
            return true;
        }
        else
        {
            //we need to reparameterise the molecule
            try
            {
                params = forcefield.getParameters(FFMoleculeBase::molecule(),
                                                  map);
            }
            catch(...)
            {
                FFMoleculeBase::restore(old_mol);
                throw;
            }
            
            return true;
        }
    }
    else
        return false;
}

/** Return the differences between this molecule and 'newmol'. This returns
    the parts of this molecule that are different to 'newmol'. 
    
    \throw SireError::incompatible_error
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecule<PTNL> FFMolecule<PTNL>::getDifferences(const FFMolecule<PTNL> &newmol) const
{
    if (this->operator==(newmol))
        return FFMolecule<PTNL>();
    
    else if (this->isEmpty() or newmol.isEmpty())
        return *this;

    //ensure that they are compatible
    SireFF::detail::assertCompatible(*this, newmol);
    
    //how many CutGroups are there in this view?
    int ncgroups = FFMoleculeBase::molecule().selection().nSelectedCutGroups();

    //if there is only one group, then it has either changed or not
    if (this->params.changedAllGroups(newmol.params))
        return *this;

    //there is more than one CutGroup in this 
    //selection - get the groups that have changed
    QSet<quint32> changed_groups = this->params.getChangedGroups(newmol.params);

    if (changed_groups.count() == ncgroups)
        //all of the CutGroups have changed
        return *this;
    
    else if (changed_groups.isEmpty())
        //nothing has changed
        return FFMolecule<PTNL>();
    
    else
    {
        //only some of the CutGroups have changed - create a selection
        //of the CutGroups that changed - remember that we need to convert
        //the index into the parameters array into the CGIdx...
        AtomSelection changed_atoms = FFMoleculeBase::molecule().selection();

        QList<CGIdx> changed_cgroups;
        
        if (changed_atoms.selectedAllCutGroups())
        {
            //the indicies correspond exactly with the CGIdxs
            foreach (quint32 changed_group, changed_groups)
            {
                changed_cgroups.append( CGIdx(changed_group) );
            }
        }
        else
        {
            QList<CGIdx> selected_cgroups = changed_atoms.selectedCutGroups();
        
            foreach (quint32 changed_group, changed_groups)
            {
                //get the CGIdx of this group
                changed_cgroups.append( selected_cgroups.at(changed_group) );
            }
        }
        
        //now mask the selection so that only atoms in the changed
        //CutGroups are selected
        changed_atoms.mask(changed_cgroups);
        
        return FFMolecule<PTNL>( PartialMolecule( FFMoleculeBase::molecule().data(),
                                                  changed_atoms ),

                                 params.applyMask(changed_groups) );
    }
}

///////
/////// Implementation of FFMolecules
///////

/** Null constructor */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecules<PTNL>::FFMolecules() : FFMoleculesBase()
{}

/** Construct by converting a MoleculeGroup

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecules<PTNL>::FFMolecules(const MoleculeGroup &molgroup, PTNL &forcefield,
                               const PropertyMap &map)
                  : FFMoleculesBase(molgroup, map)
{
    mols_by_idx = ChunkedVector<Molecule>(molgroup.nMolecules());

    int i = 0;

    for (MoleculeGroup::const_iterator it = molgroup.constBegin();
         it != molgroup.constEnd();
         ++it)
    {
        mols_by_idx[i] = Molecule(*it, forcefield, map);
        ++i;
    }
}

/** Copy constructor */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecules<PTNL>::FFMolecules(const FFMolecules<PTNL> &other)
                  : FFMoleculesBase(other),
                    mols_by_idx(other.mols_by_idx)
{}

/** Destructor */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecules<PTNL>::~FFMolecules()
{}

/** Copy assignment operator */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
FFMolecules<PTNL>& FFMolecules<PTNL>::operator=(const FFMolecules<PTNL> &other)
{
    if (this != &other)
    {
        FFMoleculesBase::operator=(other);
        mols_by_idx = other.mols_by_idx;
        parameter_names = other.parameter_names;
    }
    
    return *this;
}

/** Comparison operator */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecules<PTNL>::operator==(const FFMolecules<PTNL> &other) const
{
    return this == &other or
           (mols_by_idx == other.mols_by_idx and
            parameter_names == other.parameter_names);
}

/** Comparison operator */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecules<PTNL>::operator!=(const FFMolecules<PTNL> &other) const
{
    return not this->operator==(other);
}

/** Completely remove all of the molecules */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
void FFMolecules<PTNL>::clear()
{
    mols_by_idx.clear();
    FFMoleculesBase::clear();
}

/** Change the molecule 'new_molecule' (which is in the forcefield 'forcefield').
    This does nothing if this molecule is not in this group. If 
    'record_changes' is true then this returns a ChangedMolecule that
    records the change from the current version of the molecule 
    to 'new_molecule'. If there is no change, or record_changes is false,
    then a null (isEmpty()) ChangedMolecule is returned

    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
    \throw SireBase::missing_property
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
typename FFMolecules<PTNL>::ChangedMolecule 
FFMolecules<PTNL>::change(const SireMol::Molecule &new_molecule,
                          PTNL &forcefield, bool record_changes)
{
    QHash<MolNum,quint32>::const_iterator it = FFMoleculesBase::indexesByMolNum()
                                                    .constFind(new_molecule.number());

    if (it == FFMoleculesBase::indexesByMolNum().constEnd())
        //this molecule is not in this group
        return ChangedMolecule();

    if (record_changes)
    {
        //get the current copy of the molecule
        typename FFMolecules<PTNL>::Molecule old_molecule = mols_by_idx[*it];
        
        //now make the change
        if (mols_by_idx[*it].change(new_molecule, forcefield,
                                    parameter_names.constData()[*it]))
        {
            return ChangedMolecule(old_molecule, mols_by_idx[*it]);
        }
        else
        {
            //there was no change
            return ChangedMolecule();
        }
    }
    else
    {
        mols_by_idx[*it].change(new_molecule, forcefield,
                                parameter_names.constData()[*it]);

        return ChangedMolecule();
    }
}

/** Add the molecule 'molecule' to this group, getting the parameters
    for the molecule using the parameter names listed in 'paramnames'.
    
    If this molecule already exists in this group then this will
    add any parts of 'molecule' that don't already exist, but this
    will *not* update the existing molecule so that it is at the 
    same version. If you want to add to the same version, then
    you must call change(molecule, ffield) then add(molecule, ffield).
    
    Note also that you *cannot* add this molecule if the existing
    version has a different molecule layout ID.
    
    Note also that you *cannot* change the source of the parameters
    if this molecule exists already in the forcefield. To change
    parameter sources for a molecule that is already in the forcefield
    you must remove the molecule and then add it from scratch.

    If 'record_changes' is true then this returns a ChangedMolecule
    that records the change. If there is no change, or if 'record_changes'
    is false, then a null ChangedMolecule is returned
    
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
    \throw SireBase::missing_property
*/                                                        
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
typename FFMolecules<PTNL>::ChangedMolecule 
FFMolecules<PTNL>::add(const PartialMolecule &molecule,
                       const PropertyMap &map,
                       PTNL &forcefield, bool record_changes)
{
    //does this molecule exist already?
    if (FFMoleculesBase::contains(molecule.number()))
    {
        //difficult case - the existing molecule may be incompatible.
        quint32 i = FFMoleculesBase::indexOf(molecule.number());
        
        //ensure that the parameter source names are the same
        if (map != parameter_names.constData()[i])
            SireFF::detail::throwFFMoleculesIncompatibleParameterNames(
                    parameter_names.constData()[i], map);
        
        if (record_changes)
        {
            //get the existing copy of the molecule
            typename FFMolecules<PTNL>::Molecule old_mol = mols_by_idx[i];
            
            //try to make the change
            if (mols_by_idx[i].add(molecule.selection(), forcefield, map))
            {
                //there was a change!
                return ChangedMolecule(old_mol, mols_by_idx[i]);
            }
            else
                //there was no change
                return ChangedMolecule();
        }
        else
        {
            //add the molecule
            mols_by_idx[i].add(molecule.selection(), forcefield, map);
            
            return ChangedMolecule();
        }
    }
    else
    {
        //we are adding this molecule from scratch. First create
        //a forcefield specialised version of this molecule
        typename FFMolecules<PTNL>::Molecule ffmol(molecule, forcefield, map);
        
        //that was successful! Add the molecule into the index
        FFMoleculesBase::_pvt_add(molecule.number());
        
        mols_by_idx.append(ffmol);
        parameter_names.append(map);
        
        if (record_changes)
            return ChangedMolecule( typename FFMolecules<PTNL>::Molecule(),
                                    ffmol );
        else
            return ChangedMolecule();
    }
}
            
/** Remove the molecule 'molecule' from this forcefield. This only
    removes the selected atoms in 'molecule', and this molecule must
    have the same molecule layout version. Note that this does not
    update the version of the molecule in this set to be the 
    same as 'molecule'
    
    If 'record_changes' is true then this records the change
    caused by the removal of the atoms and returns it in the
    returned ChangedMolecule. If there is no change, or if
    'record_changes' is false, then a null ChangedMolecule
    is returned.
    
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
    \throw SireBase::missing_property
*/        
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
typename FFMolecules<PTNL>::ChangedMolecule 
FFMolecules<PTNL>::remove(const PartialMolecule &molecule,
                          PTNL &forcefield, bool record_changes)
{
    if (not this->contains(molecule.number()))
        return ChangedMolecule();

    quint32 i = this->indexOf(molecule.number());
    
    if (record_changes)
    {
        typename FFMolecules<PTNL>::Molecule old_mol = mols_by_idx[i];
        
        if (mols_by_idx[i].remove(molecule.selection(), forcefield,
                                  parameter_names.constData()[i]))
        {
            //some atoms were removed - were all of them?
            if (mols_by_idx[i].isEmpty())
            {
                //yes - the molecule has been completely removed
                FFMoleculesBase::_pvt_remove(molecule.number());
                
                mols_by_idx.remove(i);
                parameter_names.remove(i);
                
                return ChangedMolecule(old_mol, 
                                       typename FFMolecules<PTNL>::Molecule());
            }
            
            return ChangedMolecule(old_mol, mols_by_idx[i]);
        }
        else
            //nothing was removed
            return ChangedMolecule();
    }
    else
    {
        if (mols_by_idx[i].remove(molecule.selection(), forcefield,
                                  parameter_names.constData()[i]))
        {
            //has the whole molecule been removed?
            if (mols_by_idx[i].isEmpty())
            {
                //yes!
                FFMolecules::_pvt_remove(molecule.number());
                mols_by_idx.remove(i);
                parameter_names.remove(i);
            }
        }
        
        return ChangedMolecule();
    }
}

/** Return whether or not the use of the property map 'map' would
    change the properties used to get the parameters of the molecule
    with number 'molnum'
    
    \throw SireMol::missing_molecule
*/
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
bool FFMolecules<PTNL>::wouldChangeProperties(MolNum molnum,
                                              const PropertyMap &map) const
{
    return parameter_names.at(this->indexOf(molnum)) != map;
}

/** Return the array of all of the molecules in this group arranged
    in the order they appear in this group */
template<class PTNL>
SIRE_OUTOFLINE_TEMPLATE
const ChunkedVector<typename FFMolecules<PTNL>::Molecule>&
FFMolecules<PTNL>::moleculesByIndex() const
{
    return mols_by_idx;
}

///////
/////// Implementation of ChangedMolecule
///////

/** Null constructor */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>::ChangedMolecule()
{}

/** Construct from just the old molecule. This represents no change */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>::ChangedMolecule(const FFMOL &molecule)
                       : old_molecule(molecule), new_molecule(molecule),
                         old_parts(molecule), new_parts(molecule)
{}

/** Construct the change from 'oldmol' to 'newmol'. Both 'oldmol'
    and 'newmol' must be views of the same molecule (i.e. their
    molecule numbers must be identical!)

    \throw SireError::incompatible_error
*/
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>::ChangedMolecule(const FFMOL &oldmol, const FFMOL &newmol)
                       : old_molecule(oldmol),
                         new_molecule(newmol)
{
    if (oldmol.isEmpty())
    {
        //this represents creating newmol
        new_parts = new_molecule;
    }
    else if (newmol.isEmpty())
    {
        //this represent destroying oldmol
        old_parts = old_molecule;
    }
    else
    {
        //we are changing from oldmol to newmol
        SireFF::detail::assertCompatible(old_molecule, new_molecule);
        
        //now work out the differences
        if (oldmol != newmol)
        {
            old_parts = old_molecule.getDifferences(new_molecule);
            new_parts = new_molecule.getDifferences(old_molecule);
        }
    }
}

/** Copy constructor */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>::ChangedMolecule(const ChangedMolecule<FFMOL> &other)
                       : old_molecule(other.old_molecule),
                         new_molecule(other.new_molecule),
                         old_parts(other.old_parts),
                         new_parts(other.new_parts)
{}

/** Destructor */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>::~ChangedMolecule()
{}

/** Copy assignment operator */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>& ChangedMolecule<FFMOL>::operator=(
                                            const ChangedMolecule<FFMOL> &other)
{
    if (this != &other)
    {
        old_molecule = other.old_molecule;
        new_molecule = other.new_molecule;
        old_parts = other.old_parts;
        new_parts = other.new_parts;
    }
    
    return *this;
}

/** Comparison operator */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
bool ChangedMolecule<FFMOL>::operator==(const ChangedMolecule<FFMOL> &other) const
{
    return this == &other or
           (old_molecule == other.old_molecule and 
            new_molecule == other.new_molecule and
            old_parts == other.old_parts and
            new_parts == other.new_parts);
}

/** Comparison operator */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
bool ChangedMolecule<FFMOL>::operator!=(const ChangedMolecule<FFMOL> &other) const
{
    return not this->operator==(other);
}

/** Return the number of this molecule */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
MolNum ChangedMolecule<FFMOL>::number() const
{
    if (old_molecule.isEmpty())
        return new_molecule.number();
    else
        return old_molecule.number();
}

/** Return whether this is an empty change (contains nothing!) */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
bool ChangedMolecule<FFMOL>::isEmpty() const
{
    return old_molecule.isEmpty() and new_molecule.isEmpty();
}

/** Return whether or not the entire molecule has changed */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
bool ChangedMolecule<FFMOL>::changedAll() const
{
    return new_molecule == new_parts and 
           old_molecule == old_parts;
}

/** Return whether nothing has changed */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
bool ChangedMolecule<FFMOL>::nothingChanged() const
{
    return new_molecule == old_molecule;
}

/** Return the complete old version of the molecule
    (this may be empty if this represents the creation of a molecule) */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
const FFMOL& ChangedMolecule<FFMOL>::oldMolecule() const
{
    return old_molecule;
}

/** Return the complete new version of the molecule
    (this may be empty if this represents the removal of a molecule) */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
const FFMOL& ChangedMolecule<FFMOL>::newMolecule() const
{
    return new_molecule;
}

/** Return the parts of the old molecule that have changed */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
const FFMOL& ChangedMolecule<FFMOL>::oldParts() const
{
    return old_parts;
}

/** Return the parts of the new molecule that have changed */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
const FFMOL& ChangedMolecule<FFMOL>::newParts() const
{
    return new_parts;
}

/** Record that this molecule has in fact been removed */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>& ChangedMolecule<FFMOL>::removed()
{
    new_molecule = FFMOL();
    new_parts = FFMOL();
    old_parts = old_molecule;
    
    return *this;
}

/** Use this function to construct a ChangedMolecule that represents
    the removal of 'mol' */
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL> ChangedMolecule<FFMOL>::removeMolecule(const FFMOL &mol)
{
    ChangedMolecule<FFMOL> change;
    
    change.old_molecule = mol;
    change.old_parts = mol;
    
    return change;
}


/** Further change this change so that it equals the change
    from old_molecule to new_new_molecule
    
    \throw SireError::incompatible_error
*/
template<class FFMOL>
SIRE_OUTOFLINE_TEMPLATE
ChangedMolecule<FFMOL>& ChangedMolecule<FFMOL>::change(
                                            const FFMOL &new_new_molecule)
{
    if (new_new_molecule.isEmpty())
    {
        //we are now removing the molecule
        new_molecule = new_new_molecule;
        new_parts = FFMOL();
        old_parts = old_molecule;
    }
    else if (old_molecule.isEmpty())
    {
        //we are changing an added molecule - ensure that it is compatible
        if (not new_molecule.isEmpty())
            SireFF::detail::assertCompatible(new_molecule, new_new_molecule);

        new_molecule = new_new_molecule;
        new_parts = new_new_molecule;
    }
    else
    {
        SireFF::detail::assertCompatible(old_molecule, new_molecule);
        
        new_molecule = new_new_molecule;
        
        if (new_molecule != old_molecule)
        {
            old_parts = old_molecule.getDifferences(new_molecule);
            new_parts = new_molecule.getDifferences(old_molecule);
        }
        else
        {
            //this represents no change
            new_molecule = old_molecule;
            
            old_parts = FFMOL();
            new_parts = FFMOL();
        }
    }
    
    return *this;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

}

/** Serialise to a binary datastream */
template<class PTNL>
QDataStream& operator<<(QDataStream &ds, 
                        const SireFF::detail::FFMolecule<PTNL> &ffmol)
{
    SireStream::SharedDataStream sds(ds);

    sds << static_cast<const SireFF::detail::FFMoleculeBase&>(ffmol)
        << ffmol.params;
        
    return ds;
}

/** Extract from a binary datastream */
template<class PTNL>
QDataStream& operator>>(QDataStream &ds,
                        SireFF::detail::FFMolecule<PTNL> &ffmol)
{
    SireStream::SharedDataStream sds(ds);

    sds >> static_cast<SireFF::detail::FFMoleculeBase&>(ffmol)
        >> ffmol.params;
    
    return ds;
}

/** Serialise to a binary datastream */
template<class PTNL>
QDataStream& operator<<(QDataStream &ds, 
                        const SireFF::detail::FFMolecules<PTNL> &ffmols)
{
    SireStream::SharedDataStream sds(ds);

    sds << static_cast<const SireFF::detail::FFMoleculesBase&>(ffmols)
        << ffmols.mols_by_idx;
    
    return ds;
}

/** Extract from a binary datastream */
template<class PTNL>
QDataStream& operator>>(QDataStream &ds,
                        SireFF::detail::FFMolecules<PTNL> &ffmols)
{
    SireStream::SharedDataStream sds(ds);

    sds >> static_cast<SireFF::detail::FFMoleculesBase&>(ffmols)
        >> ffmols.mols_by_idx;
    
    return ds;
}

/** Serialise to a binary datastream */
template<class FFMOL>
QDataStream& operator<<(QDataStream &ds, 
                        const SireFF::detail::ChangedMolecule<FFMOL> &changedmol)
{
    SireStream::SharedDataStream sds(ds);

    sds << changedmol.old_molecule << changedmol.new_molecule
        << changedmol.old_parts << changedmol.new_parts;
    
    return ds;
}

/** Extract from a binary datastream */
template<class FFMOL>
QDataStream& operator>>(QDataStream &ds,
                        SireFF::detail::ChangedMolecule<FFMOL> &changedmol)
{
    SireStream::SharedDataStream sds(ds);

    sds >> changedmol.old_molecule >> changedmol.new_molecule
        >> changedmol.old_parts >> changedmol.new_parts;
    
    return ds;
}

SIRE_END_HEADER

#endif
