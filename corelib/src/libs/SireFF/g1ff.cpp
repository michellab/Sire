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

#include "g1ff.h"

#include "SireMol/mgnum.h"
#include "SireMol/mgname.h"
#include "SireMol/molnum.h"
#include "SireMol/molname.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireMol;
using namespace SireStream;

static const RegisterMetaType<G1FF> r_g1ff( MAGIC_ONLY, "SireFF::G1FF" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const G1FF &g1ff)
{
    writeHeader(ds, r_g1ff, 1);
    
    SharedDataStream sds(ds);

    sds << g1ff.molgroup << g1ff.allow_overlap_of_atoms
        << static_cast<const FF&>(g1ff);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      G1FF &g1ff)
{
    VersionID v = readHeader(ds, r_g1ff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> g1ff.molgroup >> g1ff.allow_overlap_of_atoms
            >> static_cast<FF&>(g1ff);
    }
    else
        throw version_error(v, "1", r_g1ff, CODELOC);
        
    return ds;
}

/** Constructor */
G1FF::G1FF(bool allow_overlap) 
     : FF(), allow_overlap_of_atoms(allow_overlap)
{
    molgroup = FFMolGroupPvt(this->name(), 0, this);
    molgroup.setParent(this);
    MolGroupsBase::addToIndex(molgroup);
}

/** Copy constructor */
G1FF::G1FF(const G1FF &other) 
     : FF(other), molgroup(other.molgroup),
       allow_overlap_of_atoms(other.allow_overlap_of_atoms)
{
    molgroup.setParent(this);
    molgroup.setIndex(0);
}

/** Copy assignment operator */
G1FF& G1FF::operator=(const G1FF &other)
{
    molgroup = other.molgroup;
    molgroup.setParent(this);
    molgroup.setIndex(0);
    
    allow_overlap_of_atoms = other.allow_overlap_of_atoms;
    FF::operator=(other);
    
    return *this;
}

/** Destructor */
G1FF::~G1FF()
{}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& G1FF::at(MGNum mgnum) const
{
    if (molgroup.number() != mgnum)
        throw SireMol::missing_group( QObject::tr(
            "The forcefield %1 (%2) does not contain molecule group "
            "with number %3. The only group it contains has number "
            "%4.")
                .arg(this->name()).arg(this->what())
                .arg(mgnum).arg(molgroup.number()), CODELOC );
                
    return molgroup;
}

/** Reindex the collection of molecule groups */
void G1FF::reindex()
{
    MolGroupsBase::clearIndex();
    MolGroupsBase::addToIndex( molgroup );
}

/** Update the name of the molecule group in this forcefield so that
    it matches the name of the forcefield */
void G1FF::_pvt_updateName()
{
    if (molgroup.name().value() == this->name().value())
        return;

    molgroup.setName( this->name() );
    molgroup.setNewNumber();
    
    this->reindex();
}

/** Assert that this forcefield contains the group with number 'mgnum'

    \throw SireMol::missing_group
*/
void G1FF::assertContains(MGNum mgnum) const
{
    if (molgroup.number() != mgnum)
        throw SireMol::missing_group( QObject::tr(
            "The forcefield %1 does not contain a group with "
            "number %2. The only group it contains has number %4.")
                .arg(this->name())
                .arg(mgnum).arg(molgroup.number()), CODELOC );
}

/** Return a reference to the group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& G1FF::getGroup(MGNum mgnum) const
{
    G1FF::assertContains(mgnum);
    return molgroup;
}

/** Return const pointers to the groups with number 'mgnums'

    \throw SireMol::missing_group
*/
void G1FF::getGroups(const QList<MGNum> &mgnums,
                     QVarLengthArray<const MoleculeGroup*,10> &groups) const
{
    groups.clear();
    
    foreach (MGNum mgnum, mgnums)
    {
        G1FF::assertContains(mgnum);
        groups.append( &molgroup );
    }
}

/** Return pointers to all of the groups in this forcefield */
QHash<MGNum,const MoleculeGroup*> G1FF::getGroups() const
{
    QHash<MGNum,const MoleculeGroup*> groups;
    groups.reserve(1);
    groups.insert( molgroup.number(), &molgroup );
    
    return groups;
}

/** Assert that i refers to a valid group

    \throw SireError::program_bug
*/
void G1FF::assertValidGroup(quint32 i) const
{
    if (i != 0)
        throw SireError::program_bug( QObject::tr(
            "G1FF should only ever use group index 0 - there is a bug "
            "somewhere as we've just been passed group index %1.")
                .arg(i), CODELOC );
}

/** Set the name of the group in this forcefield */
void G1FF::group_setName(quint32 i, const QString &new_name)
{
    assertValidGroup(i);
    molgroup.setName(new_name);
    this->reindex();
}

/** Assert that there is no overlap between the atoms in 
    'molview' and the atoms in 'group'
    
    \throw SireMol::duplicate_atom
*/
void G1FF::assertNoOverlap(const MoleculeGroup &group,
                           const MoleculeView &molview) const
{
    if (group.intersects(molview))
        throw SireMol::duplicate_atom( QObject::tr(
            "Some of the atoms in the view of the molecule %1 (%2) "
            "are already present in the forcefield group %3 (%4).")
                .arg(molview.data().name()).arg(molview.data().number())
                .arg(group.name()).arg(group.number()), CODELOC );
}

/** Tell the derived forcefield that the following views have just
    been added - use the supplied map to get the parameters */
void G1FF::_pvt_added(const ViewsOfMol &molviews, const PropertyMap &map)
{
    this->_pvt_added(molviews.all(), map);
}

/** Tell the derived forcefield that the following views have
    just been removed */
void G1FF::_pvt_removed(const ViewsOfMol &molviews)
{
    this->_pvt_removed(molviews.all());
}

/** Record that all copies of the view in 'mol' have been removed */
void G1FF::_pvt_removedAll(const PartialMolecule &mol)
{
    this->_pvt_removed(mol);
}

/** Record that all copies of the views in 'mol' have been removed */
void G1FF::_pvt_removedAll(const ViewsOfMol &mol)
{
    this->_pvt_removed(mol);
}

/** Assert that there is no overlap between the atoms in 'molview'
    and any atoms that exist currently in this molecule
    
    \throw SireMol::duplicate_atom
*/
void G1FF::assertNoOverlap(const MoleculeView &molview) const
{
    this->assertNoOverlap(molgroup, molview);
}

/** Add the molecule view in 'molview' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void G1FF::group_add(quint32 i, const MoleculeView &molview,
                     const PropertyMap &map)
{
    assertValidGroup(i);

    if (molview.selection().isEmpty())
        return;

    if (not allow_overlap_of_atoms)
        this->assertNoOverlap(molview);

    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        molgroup.add(molview);
        this->_pvt_added( PartialMolecule(molview), map );
        
        FF::incrementVersion();
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
}

/** Add the views of the molecule in 'molviews' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void G1FF::group_add(quint32 i, const ViewsOfMol &molviews, 
                     const PropertyMap &map)
{
    assertValidGroup(i);

    if (molviews.isEmpty())
        return;

    if (not allow_overlap_of_atoms)
    {
        this->assertNoOverlap(molviews);
        molviews.assertNoOverlap();
    }
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        molgroup.add(molviews);
        this->_pvt_added(molviews, map);
        FF::incrementVersion();
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
}

/** Add the molecules in 'molecules' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void G1FF::group_add(quint32 i, const Molecules &molecules, 
                     const PropertyMap &map)
{
    assertValidGroup(i);

    if (molecules.isEmpty())
        return;

    if (not allow_overlap_of_atoms)
    {
        //assert that there is no overlap between the molecules to be added
        //and the molecules already existing in this forcefield
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            this->assertNoOverlap(*it);
        
            //also assert that there is no overlap within the molecule
            it->assertNoOverlap();
        }
    }
    
    //save the old state of this forcefield
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        if (not allow_overlap_of_atoms)
        {
            //add the molecules - they must be unique views of we'd
            //have thrown an exception earlier!
            QList<ViewsOfMol> added_mols = molgroup.addIfUnique(molecules);
            
            foreach (const ViewsOfMol &added_mol, added_mols)
            {
                this->_pvt_added(added_mol, map);
            }
        }
        else
        {
            //add the molecules...
            molgroup.add(molecules);
            
            //now convert the 'Molecules' into a list...
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                this->_pvt_added(*it, map);
            }
        }
        
        FF::incrementVersion();
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
}

/** Add the views of the molecules in the group 'molgroup' to this forcefield, 
    using the supplied property map to get the names of the properties that 
    contain the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void G1FF::group_add(quint32 i, const MoleculeGroup &new_group, 
                     const PropertyMap &map)
{
    G1FF::group_add(i, new_group.molecules(), map);
}

/** Add the molecule view in 'molview' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool G1FF::group_addIfUnique(quint32 i, const MoleculeView &molview, 
                             const PropertyMap &map)
{
    assertValidGroup(i);

    if (molview.selection().isEmpty())
        return false;

    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        if (molgroup.addIfUnique(molview))
        {
            if (not allow_overlap_of_atoms)
            {
                //the molecule view was added successfully
                // - we must ensure that this view did not overlap
                //   with any of the existing atoms
                this->assertNoOverlap(old_state, molview);
            }

            //now rebuild this molecule in the forcefield
            this->_pvt_added( PartialMolecule(molview), map );

            FF::incrementVersion();
            
            return true;
        }
        else
            return false;
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return false;
}
                             
/** Add the views of the molecule in 'molviews' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
ViewsOfMol G1FF::group_addIfUnique(quint32 i, const ViewsOfMol &molviews, 
                                   const PropertyMap &map)
{
    assertValidGroup(i);

    if (molviews.isEmpty())
        return ViewsOfMol();
          
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        ViewsOfMol added_views = molgroup.addIfUnique(molviews);
        
        if (not added_views.isEmpty())
        {
            if (not allow_overlap_of_atoms)
            {
                this->assertNoOverlap(old_state, added_views);
                added_views.assertNoOverlap();
            }

            this->_pvt_added(added_views, map);
            
            FF::incrementVersion();
        }
        
        return added_views;
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Add the molecules in 'molecules' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
QList<ViewsOfMol> G1FF::group_addIfUnique(quint32 i, const Molecules &molecules, 
                                          const PropertyMap &map)
{
    assertValidGroup(i);

    if (molecules.isEmpty())
        return QList<ViewsOfMol>();
        
    FFMolGroupPvt old_molgroup = molgroup;
        
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> added_mols = molgroup.addIfUnique(molecules);
        
        if (not added_mols.isEmpty())
        {
            if (not allow_overlap_of_atoms)
            {
                //assert that there is no overlap between the added molecules
                //and the existing molecules
                foreach (const ViewsOfMol &added_mol, added_mols)
                {
                    this->assertNoOverlap(old_molgroup, added_mol);
                    added_mol.assertNoOverlap();
                }
            }
            
            //now get the parameters
            foreach (const ViewsOfMol &added_mol, added_mols)
            {
                this->_pvt_added(added_mol, map);
            }
            
            FF::incrementVersion();
        }
        
        return added_mols;
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return QList<ViewsOfMol>();
}

/** Add the views of the molecules in the group 'molgroup' to this forcefield, 
    using the supplied property map to get the names of the properties that 
    contain the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
QList<ViewsOfMol> G1FF::group_addIfUnique(quint32 i, const MoleculeGroup &new_group, 
                                          const PropertyMap &map)
{
    return G1FF::group_addIfUnique(i, new_group.molecules(), map);
}

/** Remove the molecule view in 'molview' from this forcefield */
bool G1FF::group_remove(quint32 i, const MoleculeView &molview)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        if (molgroup.remove(molview))
        {
            this->_pvt_removed( PartialMolecule(molview) );
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return false;
}

/** Remove all of the views in 'molviews' from this forcefield */
ViewsOfMol G1FF::group_remove(quint32 i, const ViewsOfMol &molviews)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        ViewsOfMol removed_views = molgroup.remove(molviews);
        
        if (not removed_views.isEmpty())
        {
            this->_pvt_removed(removed_views);
            
            FF::incrementVersion();
        }
            
        return removed_views;
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Remove the molecules in 'molecules' from this forcefield */
QList<ViewsOfMol> G1FF::group_remove(quint32 i, const Molecules &molecules)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> removed_mols = molgroup.remove(molecules);
        
        if (not removed_mols.isEmpty())
        {
            foreach (const ViewsOfMol &removed_mol, removed_mols)
            {
                this->_pvt_removed(removed_mol);
            }
            
            FF::incrementVersion();
        }
            
        return removed_mols;
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return QList<ViewsOfMol>();
}

/** Remove all of the molecules in the molecule group 'molgroup' from this
    forcefield */
QList<ViewsOfMol> G1FF::group_remove(quint32 i, const MoleculeGroup &new_group)
{
    return G1FF::group_remove(i, new_group.molecules());
}

/** Remove the molecule view in 'molview' from this forcefield */
bool G1FF::group_removeAll(quint32 i, const MoleculeView &molview)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        if (molgroup.removeAll(molview))
        {
            this->_pvt_removedAll( PartialMolecule(molview) );
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return false;
}

/** Remove all of the views in 'molviews' from this forcefield */
ViewsOfMol G1FF::group_removeAll(quint32 i, const ViewsOfMol &molviews)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        ViewsOfMol removed_views = molgroup.removeAll(molviews);
        
        if (not removed_views.isEmpty())
        {
            this->_pvt_removedAll(removed_views);
            
            FF::incrementVersion();
        }
            
        return removed_views;
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Remove the molecules in 'molecules' from this forcefield */
QList<ViewsOfMol> G1FF::group_removeAll(quint32 i, const Molecules &molecules)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> removed_mols = molgroup.removeAll(molecules);
        
        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            this->_pvt_removedAll(removed_mol);
            
            FF::incrementVersion();
        }
            
        return removed_mols;
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return QList<ViewsOfMol>();
}

/** Remove all of the molecules in the molecule group 'molgroup' from this
    forcefield */
QList<ViewsOfMol> G1FF::group_removeAll(quint32 i, const MoleculeGroup &new_group)
{
    return G1FF::group_removeAll(i, new_group.molecules());
}

/** Remove the molecule with number 'molnum' from this forcefield */
ViewsOfMol G1FF::group_remove(quint32 i, MolNum molnum)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        ViewsOfMol removed_mol = molgroup.remove(molnum);
        
        if (not removed_mol.isEmpty())
        {
            this->_pvt_removed(removed_mol);
            
            FF::incrementVersion();
        }
            
        return removed_mol;
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Remove all of the molecules whose numbers are in 'molnum' from
    this forcefield */
QList<ViewsOfMol> G1FF::group_remove(quint32 i, const QSet<MolNum> &molnums)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> removed_mols = molgroup.remove(molnums);
        
        if (not removed_mols.isEmpty())
        {
            foreach (const ViewsOfMol &removed_mol, removed_mols)
            {
                this->_pvt_removed(removed_mol);
            }
            
            FF::incrementVersion();
        }
            
        return removed_mols;
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return QList<ViewsOfMol>();
}

/** Completely remove all molecules from this forcefield */
void G1FF::group_removeAll(quint32 i)
{
    molgroup.removeAll();
    this->_pvt_removedAll();
}

/** Update the molecule whose data is in 'moldata' to use this
    version of the molecule data in this forcefield */
bool G1FF::group_update(quint32 i, const MoleculeData &moldata, bool auto_commit)
{
    assertValidGroup(i);

    if (molgroup.update(moldata, auto_commit))
    {
        this->_pvt_changed( Molecule(moldata), auto_commit );
        
        FF::incrementVersion();
        
        return true;
    }
    else
        return false;
}

/** Update this forcefield so that it uses the same version of the
    molecules as in 'molecules' */
QList<Molecule> G1FF::group_update(quint32 i, const Molecules &molecules, bool auto_commit)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup;
    
    try
    {
        QList<Molecule> updated_mols = molgroup.update(molecules, auto_commit);
        
        if (not updated_mols.isEmpty())
        {
            this->_pvt_changed(updated_mols, auto_commit);
            
            FF::incrementVersion();
        }
            
        return updated_mols;
    }
    catch(...)
    {
        molgroup = old_state;
        throw;
    }
    
    return QList<Molecule>();
}

/** Update the molecule group in this forcefield so that it has
    the same molecule versions as in 'new_group' */
QList<Molecule> G1FF::group_update(quint32 i, const MoleculeGroup &new_group, bool auto_commit)
{
    return G1FF::group_update(i, new_group.molecules(), auto_commit);
}

/** Set the contents of this forcefield so that it only contains the 
    view of the molecule in 'molview' (using the supplied property
    map to find the properties that contain the required parameters
    for this forcefield) */
bool G1FF::group_setContents(quint32 i, const MoleculeView &molview, 
                             const PropertyMap &map)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup.setContents(molview);
        
        if (changed or this->_pvt_wouldChangeProperties(molview.data().number(),map))
        {
            this->_pvt_removedAll();
            this->_pvt_added( PartialMolecule(molview), map );
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return false;
}

/** Set the contents of this forcefield so that it contains just
    the views of the molecule in 'molviews' */
bool G1FF::group_setContents(quint32 i, const ViewsOfMol &molviews, 
                             const PropertyMap &map)
{
    assertValidGroup(i);
    
    if (not allow_overlap_of_atoms)
        molviews.assertNoOverlap();
        
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup.setContents(molviews);
        
        if (changed or this->_pvt_wouldChangeProperties(molviews.number(), map))
        {
            this->_pvt_removedAll();
            this->_pvt_added(molviews, map);
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return false;
}

/** Set the contents of this forcefield so that it only contains 'molecules' */
bool G1FF::group_setContents(quint32 i, const Molecules &molecules, 
                             const PropertyMap &map)
{
    assertValidGroup(i);
    
    if (not allow_overlap_of_atoms)
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            it->assertNoOverlap();
        }
    }
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup.setContents(molecules);
        
        if (not changed)
        {
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (this->_pvt_wouldChangeProperties(it->number(), map))
                {
                    changed = true;
                    break;
                }
            }
        }
        
        if (changed)
        {
            this->_pvt_removedAll();
            
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                this->_pvt_added(*it, map);
            }
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return false;
}

/** Set the contents of this forcefield so that it contains only
    the molecules in 'molgroup' */
bool G1FF::group_setContents(quint32 i, const MoleculeGroup &new_group, 
                             const PropertyMap &map)
{
    assertValidGroup(i);
    
    if (not allow_overlap_of_atoms)
    {
        for (MoleculeGroup::const_iterator it = new_group.constBegin();
             it != new_group.constEnd();
             ++it)
        {
            it->assertNoOverlap();
        }
    }
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup.setContents(new_group);
        
        if (not changed)
        {
            for (MoleculeGroup::const_iterator it = new_group.constBegin();
                 it != new_group.constEnd();
                 ++it)
            {
                if (this->_pvt_wouldChangeProperties(it->number(), map))
                {
                    changed = true;
                    break;
                }
            }
        }
        
        if (changed)
        {
            this->_pvt_removedAll();
            
            for (MoleculeGroup::const_iterator it = new_group.constBegin();
                 it != new_group.constEnd();
                 ++it)
            {
                this->_pvt_added(*it, map);
            }
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
    
    return false;
}

void G1FF::add(const MoleculeView &molview, const PropertyMap &map)
{
    FF::add(molview, MGIdx(0), map);
}

void G1FF::add(const ViewsOfMol &molviews, const PropertyMap &map)
{
    FF::add(molviews, MGIdx(0), map);
}

void G1FF::add(const Molecules &molecules, const PropertyMap &map)
{
    FF::add(molecules, MGIdx(0), map);
}

void G1FF::add(const MoleculeGroup &group, const PropertyMap &map)
{
    FF::add(group, MGIdx(0), map);
}

void G1FF::addIfUnique(const MoleculeView &molview, const PropertyMap &map)
{
    FF::addIfUnique(molview, MGIdx(0), map);
}

void G1FF::addIfUnique(const ViewsOfMol &molviews, const PropertyMap &map)
{
    FF::addIfUnique(molviews, MGIdx(0), map);
}

void G1FF::addIfUnique(const Molecules &molecules, const PropertyMap &map)
{
    FF::addIfUnique(molecules, MGIdx(0), map);
}

void G1FF::addIfUnique(const MoleculeGroup &group, const PropertyMap &map)
{
    FF::addIfUnique(group, MGIdx(0), map);
}

void G1FF::add(const MoleculeView &molview)
{
    FF::add(molview, MGIdx(0));
}

void G1FF::add(const ViewsOfMol &molviews)
{
    FF::add(molviews, MGIdx(0));
}

void G1FF::add(const Molecules &molecules)
{
    FF::add(molecules, MGIdx(0));
}

void G1FF::add(const MoleculeGroup &group)
{
    FF::add(group, MGIdx(0));
}

void G1FF::addIfUnique(const MoleculeView &molview)
{
    FF::addIfUnique(molview, MGIdx(0));
}

void G1FF::addIfUnique(const ViewsOfMol &molviews)
{
    FF::addIfUnique(molviews, MGIdx(0));
}

void G1FF::addIfUnique(const Molecules &molecules)
{
    FF::addIfUnique(molecules, MGIdx(0));
}

void G1FF::addIfUnique(const MoleculeGroup &group)
{
    FF::addIfUnique(group, MGIdx(0));
}

bool G1FF::removeAll()
{
    return FF::removeAll(MGIdx(0));
}

bool G1FF::remove(const MoleculeView &molview)
{
    return FF::remove(molview, MGIdx(0));
}

bool G1FF::remove(const ViewsOfMol &molviews)
{
    return FF::remove(molviews, MGIdx(0));
}

bool G1FF::remove(const Molecules &molecules)
{
    return FF::remove(molecules, MGIdx(0));
}

bool G1FF::remove(const MoleculeGroup &group)
{
    return FF::remove(group, MGIdx(0));
}

bool G1FF::removeAll(const MoleculeView &molview)
{
    return FF::removeAll(molview, MGIdx(0));
}

bool G1FF::removeAll(const ViewsOfMol &molviews)
{
    return FF::removeAll(molviews, MGIdx(0));
}

bool G1FF::removeAll(const Molecules &molecules)
{
    return FF::removeAll(molecules, MGIdx(0));
}

bool G1FF::removeAll(const MoleculeGroup &group)
{
    return FF::removeAll(group, MGIdx(0));
}

bool G1FF::remove(MolNum molnum)
{
    return FF::remove(molnum, MGIdx(0));
}

bool G1FF::remove(const QSet<MolNum> &molnums)
{
    return FF::remove(molnums, MGIdx(0));
}

void G1FF::setContents(const MoleculeView &molview)
{
    FF::setContents(MGIdx(0), molview);
}

void G1FF::setContents(const ViewsOfMol &molviews)
{
    FF::setContents(MGIdx(0), molviews);
}

void G1FF::setContents(const Molecules &molecules)
{
    FF::setContents(MGIdx(0), molecules);
}

void G1FF::setContents(const MoleculeGroup &group)
{
    FF::setContents(MGIdx(0), group);
}

void G1FF::setContents(const MoleculeView &molview, const PropertyMap &map)
{
    FF::setContents(MGIdx(0), molview, map);
}

void G1FF::setContents(const ViewsOfMol &molviews, const PropertyMap &map)
{
    FF::setContents(MGIdx(0), molviews, map);
}

void G1FF::setContents(const Molecules &molecules, const PropertyMap &map)
{
    FF::setContents(MGIdx(0), molecules, map);
}

void G1FF::setContents(const MoleculeGroup &group, const PropertyMap &map)
{
    FF::setContents(MGIdx(0), group, map);
}

/** Return whether or not this forcefield is using temporary workspace that needs 
    to be accepted */
bool G1FF::needsAccepting() const
{
    return molgroup.needsAccepting();
}

/** Tell the forcefield that the last move was accepted. This tells the
    forcefield to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void G1FF::accept()
{
    if (molgroup.needsAccepting())
        molgroup.accept();
}
