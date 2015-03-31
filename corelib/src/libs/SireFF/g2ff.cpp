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

#include "g2ff.h"

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

namespace SireFF
{
namespace detail
{

void SIREFF_EXPORT throwIntra2B2GFFIncompatibeScaleFactorsError(
                const QString &molname, MolNum molnum, quint64 v0, quint64 v1)
{
    throw SireError::incompatible_error( QObject::tr(
        "Cannot add or update the molecule %1 (%2) as the intramolecular "
        "non-bonded scaling factors are different in the added or "
        "changed molecule (version %3) to those in the molecule that "
        "is already present in this forcefield (version %4).")
            .arg(molname).arg(molnum).arg(v0).arg(v1), CODELOC );
}

} // end of namespace detail
} // end of namespace SireFF

static const RegisterMetaType<G2FF> r_g2ff( MAGIC_ONLY, "SireFF::G2FF" );

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds,
                                      const G2FF &g2ff)
{
    writeHeader(ds, r_g2ff, 1);
    
    SharedDataStream sds(ds);
    
    sds << g2ff.molgroup[0] << g2ff.molgroup[1] 
        << g2ff.allow_overlap_of_atoms
        << static_cast<const FF&>(g2ff);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds,
                                      G2FF &g2ff)
{
    VersionID v = readHeader(ds, r_g2ff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> g2ff.molgroup[0] >> g2ff.molgroup[1]
            >> g2ff.allow_overlap_of_atoms
            >> static_cast<FF&>(g2ff);
    }
    else
        throw version_error(v, "1", r_g2ff, CODELOC);
        
    return ds;
}

/** Constructor */
G2FF::G2FF(bool allow_overlap) : FF(), allow_overlap_of_atoms(allow_overlap)
{
    molgroup[0] = FFMolGroupPvt( QString("%1_A").arg(this->name()), 0, this );
    molgroup[1] = FFMolGroupPvt( QString("%1_B").arg(this->name()), 1, this );
    
    molgroup[0].setParent(this);
    molgroup[1].setParent(this);
    
    MolGroupsBase::addToIndex(molgroup[0]);
    MolGroupsBase::addToIndex(molgroup[1]);
}

/** Reindex the two groups */
void G2FF::reindex()
{
    MolGroupsBase::clearIndex();
    MolGroupsBase::addToIndex( molgroup[0] );
    MolGroupsBase::addToIndex( molgroup[1] );
}

/** Update the name of the molecule group in this forcefield so that
    it matches the name of the forcefield */
void G2FF::_pvt_updateName()
{
    QString name_a = QString("%1_A").arg(this->name());
    QString name_b = QString("%1_B").arg(this->name());

    bool need_reindex = false;

    if (molgroup[0].name().value() != name_a)
    {
        molgroup[0].setName(name_a);
        molgroup[0].setNewNumber();
        need_reindex = true;
    }
    
    if (molgroup[1].name().value() != name_b)
    {
        molgroup[1].setName(name_b);
        molgroup[1].setNewNumber();
        need_reindex = true;
    }
    
    if (need_reindex)
        this->reindex();
}

/** Copy constructor */
G2FF::G2FF(const G2FF &other) 
     : FF(other),
       allow_overlap_of_atoms(other.allow_overlap_of_atoms)
{
    molgroup[0] = other.molgroup[0];
    molgroup[1] = other.molgroup[1];
    
    molgroup[0].setParent(this);
    molgroup[1].setParent(this);
    
    molgroup[0].setIndex(0);
    molgroup[1].setIndex(1);
}

/** Destructor */
G2FF::~G2FF()
{}

/** Copy assignment operator */
G2FF& G2FF::operator=(const G2FF &other)
{
    molgroup[0] = other.molgroup[0];
    molgroup[1] = other.molgroup[1];
    
    molgroup[0].setParent(this);
    molgroup[1].setParent(this);
    molgroup[0].setIndex(0);
    molgroup[1].setIndex(1);
    
    allow_overlap_of_atoms = other.allow_overlap_of_atoms;
    FF::operator=(other);
    
    return *this;
}

/** Assert that this forcefield contains the group with number 'mgnum'

    \throw SireMol::missing_group
*/
void G2FF::assertContains(MGNum mgnum) const
{
    if (molgroup[0].number() != mgnum and molgroup[1].number() != mgnum)
        throw SireMol::missing_group( QObject::tr(
            "The forcefield %1 does not contain a group with "
            "number %2. The only groups it contains have numbers %4 and %5.")
                .arg(this->name())
                .arg(mgnum).arg(molgroup[0].number())
                .arg(molgroup[1].number()), CODELOC );
}

/** Return a reference to the group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& G2FF::getGroup(MGNum mgnum) const
{
    G2FF::assertContains(mgnum);
    
    if (molgroup[0].number() == mgnum)
        return molgroup[0];
    else
        return molgroup[1];
}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& G2FF::at(MGNum mgnum) const
{
    return this->getGroup(mgnum);
}

/** Return const pointers to the groups with number 'mgnums'

    \throw SireMol::missing_group
*/
void G2FF::getGroups(const QList<MGNum> &mgnums,
                     QVarLengthArray<const MoleculeGroup*,10> &groups) const
{
    groups.clear();
    
    foreach (MGNum mgnum, mgnums)
    {
        groups.append( &(this->at(mgnum)) );
    }
}

/** Return pointers to all of the groups in this forcefield */
QHash<MGNum,const MoleculeGroup*> G2FF::getGroups() const
{
    QHash<MGNum,const MoleculeGroup*> groups;
    groups.reserve(2);
    groups.insert( molgroup[0].number(), &(molgroup[0]) );
    groups.insert( molgroup[1].number(), &(molgroup[1]) );
    
    return groups;
}

/** Assert that i refers to a valid group

    \throw SireError::program_bug
*/
void G2FF::assertValidGroup(quint32 i) const
{
    if (i > 1)
        throw SireError::program_bug( QObject::tr(
            "G2FF should only ever use group index 0 or 1 - there is a bug "
            "somewhere as we've just been passed group index %1.")
                .arg(i), CODELOC );
}

/** Set the name of the group in this forcefield */
void G2FF::group_setName(quint32 i, const QString &new_name)
{
    assertValidGroup(i);
    molgroup[i].setName(new_name);
    this->reindex();
}

/** Assert that there is no overlap between the atoms in 
    'molview' and the atoms in 'group'
    
    \throw SireMol::duplicate_atom
*/
void G2FF::assertNoOverlap(const MoleculeGroup &group,
                           const MoleculeView &molview) const
{
    if (group.intersects(molview))
        throw SireMol::duplicate_atom( QObject::tr(
            "Some of the atoms in the view of the molecule %1 (%2) "
            "are already present in the forcefield group %3 (%4).")
                .arg(molview.data().name()).arg(molview.data().number())
                .arg(group.name()).arg(group.number()), CODELOC );
}

/** Assert that there is no overlap between the atoms in the molecules 
    in 'molecules' and the atoms in the group 'group'
    
    \throw SireMol::duplicate_atom
*/
void G2FF::assertNoOverlap(const MoleculeGroup &group,
                           const Molecules &molecules) const
{
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        this->assertNoOverlap(group, *it);
    }
}

/** Tell the derived forcefield that the following views have just
    been added to group 'groupid' - use the supplied map to get the parameters */
void G2FF::_pvt_added(quint32 groupid, const ViewsOfMol &molviews, 
                      const PropertyMap &map)
{
    this->_pvt_added(groupid, molviews.all(), map);
}

/** Tell the derived forcefield that the following views have
    just been removed */
void G2FF::_pvt_removed(quint32 groupid, const ViewsOfMol &molviews)
{
    this->_pvt_removed(groupid, molviews.all());
}

/** Record that all copies of the view in 'mol' have been removed */
void G2FF::_pvt_removedAll(quint32 groupid, const PartialMolecule &mol)
{
    this->_pvt_removed(groupid, mol);
}

/** Record that all copies of the views in 'mol' have been removed */
void G2FF::_pvt_removedAll(quint32 groupid, const ViewsOfMol &mol)
{
    this->_pvt_removed(groupid, mol);
}

/** Assert that there is no overlap between the atoms in 'molview'
    and any atoms that exist currently in this molecule
    
    \throw SireMol::duplicate_atom
*/
void G2FF::assertNoOverlap(const MoleculeView &molview) const
{
    this->assertNoOverlap(molgroup[0], molview);
    this->assertNoOverlap(molgroup[1], molview);
}

/** Add the molecule view in 'molview' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void G2FF::group_add(quint32 i, const MoleculeView &molview,
                     const PropertyMap &map)
{
    assertValidGroup(i);

    if (molview.selection().isEmpty())
        return;

    if (not allow_overlap_of_atoms)
        this->assertNoOverlap(molview);

    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        molgroup[i].add(molview);
        this->_pvt_added( i, PartialMolecule(molview), map );
        FF::incrementVersion();
    }
    catch(...)
    {
        molgroup[i] = old_state;
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
void G2FF::group_add(quint32 i, const ViewsOfMol &molviews, 
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
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        molgroup[i].add(molviews);
        this->_pvt_added(i, molviews, map);
        FF::incrementVersion();
    }
    catch(...)
    {
        molgroup[i] = old_state;
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
void G2FF::group_add(quint32 i, const Molecules &molecules, 
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
            QList<ViewsOfMol> added_mols = molgroup[i].addIfUnique(molecules);
            
            foreach (const ViewsOfMol &added_mol, added_mols)
            {
                this->_pvt_added(i, added_mol, map);
            }
        }
        else
        {
            //add the molecules...
            molgroup[i].add(molecules);
            
            //now parameterise the molecules and add them
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                this->_pvt_added(i, *it, map);
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
void G2FF::group_add(quint32 i, const MoleculeGroup &new_group, 
                     const PropertyMap &map)
{
    G2FF::group_add(i, new_group.molecules(), map);
}

/** Add the molecule view in 'molview' to this forcefield, using the 
    supplied property map to get the names of the properties that contain
    the required forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool G2FF::group_addIfUnique(quint32 i, const MoleculeView &molview, 
                             const PropertyMap &map)
{
    assertValidGroup(i);

    if (molview.selection().isEmpty())
        return false;

    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        if (molgroup[i].addIfUnique(molview))
        {
            if (not allow_overlap_of_atoms)
            {
                //the molecule view was added successfully
                // - we must ensure that this view did not overlap
                //   with any of the existing atoms
                this->assertNoOverlap(old_state, molview);
            }

            //now rebuild this molecule in the forcefield
            this->_pvt_added( i, PartialMolecule(molview), map );

            FF::incrementVersion();
            
            return true;
        }
        else
            return false;
    }
    catch(...)
    {
        molgroup[i] = old_state;
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
ViewsOfMol G2FF::group_addIfUnique(quint32 i, const ViewsOfMol &molviews, 
                                   const PropertyMap &map)
{
    assertValidGroup(i);

    if (molviews.isEmpty())
        return ViewsOfMol();
          
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        ViewsOfMol added_views = molgroup[i].addIfUnique(molviews);
        
        if (not added_views.isEmpty())
        {
            if (not allow_overlap_of_atoms)
            {
                this->assertNoOverlap(old_state, added_views);
                added_views.assertNoOverlap();
            }

            this->_pvt_added(i, added_views, map);
            
            FF::incrementVersion();
        }
        
        return added_views;
    }
    catch(...)
    {
        molgroup[i] = old_state;
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
QList<ViewsOfMol> G2FF::group_addIfUnique(quint32 i, const Molecules &molecules, 
                                          const PropertyMap &map)
{
    assertValidGroup(i);

    if (molecules.isEmpty())
        return QList<ViewsOfMol>();
        
    FFMolGroupPvt old_molgroup = molgroup[i];
        
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> added_mols = molgroup[i].addIfUnique(molecules);
        
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
                this->_pvt_added(i, added_mol, map);
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
QList<ViewsOfMol> G2FF::group_addIfUnique(quint32 i, const MoleculeGroup &new_group, 
                                          const PropertyMap &map)
{
    return G2FF::group_addIfUnique(i, new_group.molecules(), map);
}

/** Remove the molecule view in 'molview' from this forcefield */
bool G2FF::group_remove(quint32 i, const MoleculeView &molview)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        if (molgroup[i].remove(molview))
        {
            this->_pvt_removed( i, PartialMolecule(molview) );
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return false;
}

/** Remove all of the views in 'molviews' from this forcefield */
ViewsOfMol G2FF::group_remove(quint32 i, const ViewsOfMol &molviews)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        ViewsOfMol removed_views = molgroup[i].remove(molviews);
        
        if (not removed_views.isEmpty())
        {
            this->_pvt_removed(i, removed_views);
            
            FF::incrementVersion();
        }
            
        return removed_views;
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Remove the molecules in 'molecules' from this forcefield */
QList<ViewsOfMol> G2FF::group_remove(quint32 i, const Molecules &molecules)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> removed_mols = molgroup[i].remove(molecules);
        
        if (not removed_mols.isEmpty())
        {
            foreach (const ViewsOfMol &removed_mol, removed_mols)
            {
                this->_pvt_removed(i, removed_mol);
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
QList<ViewsOfMol> G2FF::group_remove(quint32 i, const MoleculeGroup &new_group)
{
    return G2FF::group_remove(i, new_group.molecules());
}

/** Remove the molecule view in 'molview' from this forcefield */
bool G2FF::group_removeAll(quint32 i, const MoleculeView &molview)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        if (molgroup[i].removeAll(molview))
        {
            this->_pvt_removedAll( i, PartialMolecule(molview) );
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return false;
}

/** Remove all of the views in 'molviews' from this forcefield */
ViewsOfMol G2FF::group_removeAll(quint32 i, const ViewsOfMol &molviews)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        ViewsOfMol removed_views = molgroup[i].removeAll(molviews);
        
        if (not removed_views.isEmpty())
        {
            this->_pvt_removedAll(i, removed_views);
            
            FF::incrementVersion();
        }
            
        return removed_views;
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Remove the molecules in 'molecules' from this forcefield */
QList<ViewsOfMol> G2FF::group_removeAll(quint32 i, const Molecules &molecules)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> removed_mols = molgroup[i].removeAll(molecules);
        
        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            this->_pvt_removedAll(i, removed_mol);
            
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
QList<ViewsOfMol> G2FF::group_removeAll(quint32 i, const MoleculeGroup &new_group)
{
    return G2FF::group_removeAll(i, new_group.molecules());
}

/** Remove the molecule with number 'molnum' from this forcefield */
ViewsOfMol G2FF::group_remove(quint32 i, MolNum molnum)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        ViewsOfMol removed_mol = molgroup[i].remove(molnum);
        
        if (not removed_mol.isEmpty())
        {
            this->_pvt_removed(i, removed_mol);
            
            FF::incrementVersion();
        }
            
        return removed_mol;
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return ViewsOfMol();
}

/** Remove all of the molecules whose numbers are in 'molnum' from
    this forcefield */
QList<ViewsOfMol> G2FF::group_remove(quint32 i, const QSet<MolNum> &molnums)
{
    assertValidGroup(i);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        QList<ViewsOfMol> removed_mols = molgroup[i].remove(molnums);
        
        if (not removed_mols.isEmpty())
        {
            foreach (const ViewsOfMol &removed_mol, removed_mols)
            {
                this->_pvt_removed(i, removed_mol);
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
void G2FF::group_removeAll(quint32 i)
{
    assertValidGroup(i);

    molgroup[i].removeAll();
    this->_pvt_removedAll(i);
}

/** Update the molecule whose data is in 'moldata' to use this
    version of the molecule data in this forcefield */
bool G2FF::group_update(quint32 i, const MoleculeData &moldata, bool auto_commit)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        if (molgroup[i].update(moldata,auto_commit))
        {
            this->_pvt_changed( i, Molecule(moldata), auto_commit );
            
            FF::incrementVersion();
            
            return true;
        }
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return false;
}

/** Update this forcefield so that it uses the same version of the
    molecules as in 'molecules' */
QList<Molecule> G2FF::group_update(quint32 i, const Molecules &molecules, bool auto_commit)
{
    assertValidGroup(i);
    
    FFMolGroupPvt old_state = molgroup[i];
    
    try
    {
        QList<Molecule> updated_mols = molgroup[i].update(molecules, auto_commit);
        
        if (not updated_mols.isEmpty())
        {
            this->_pvt_changed(i, updated_mols, auto_commit);
            
            FF::incrementVersion();
        }
            
        return updated_mols;
    }
    catch(...)
    {
        molgroup[i] = old_state;
        throw;
    }
    
    return QList<Molecule>();
}

/** Update the molecule group in this forcefield so that it has
    the same molecule versions as in 'new_group' */
QList<Molecule> G2FF::group_update(quint32 i, const MoleculeGroup &new_group, bool auto_commit)
{
    return G2FF::group_update(i, new_group.molecules(), auto_commit);
}

/** Set the contents of this forcefield so that it only contains the 
    view of the molecule in 'molview' (using the supplied property
    map to find the properties that contain the required parameters
    for this forcefield) */
bool G2FF::group_setContents(quint32 i, const MoleculeView &molview, 
                             const PropertyMap &map)
{
    assertValidGroup(i);
    
    if (not allow_overlap_of_atoms)
        this->assertNoOverlap(molgroup[i!=1], molview);
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup[i].setContents(molview);
        
        if (changed or this->_pvt_wouldChangeProperties(i, molview.data().number(),map))
        {
            this->_pvt_removedAll(i);
            this->_pvt_added( i, PartialMolecule(molview), map );
            
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
bool G2FF::group_setContents(quint32 i, const ViewsOfMol &molviews, 
                             const PropertyMap &map)
{
    assertValidGroup(i);
    
    if (not allow_overlap_of_atoms)
    {
        molviews.assertNoOverlap();
        this->assertNoOverlap(molgroup[i!=1], molviews);
    }
        
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup[i].setContents(molviews);
        
        if (changed or this->_pvt_wouldChangeProperties(i, molviews.number(), map))
        {
            this->_pvt_removedAll(i);
            this->_pvt_added(i, molviews, map);
            
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
bool G2FF::group_setContents(quint32 i, const Molecules &molecules, 
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
        
        //assert that none of these molecules overlap with the
        //molecules in the other group
        this->assertNoOverlap(molgroup[i!=1], molecules);
    }
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup[i].setContents(molecules);
        
        if (not changed)
        {
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (this->_pvt_wouldChangeProperties(i, it->number(), map))
                {
                    changed = true;
                    break;
                }
            }
        }
        
        if (changed)
        {
            this->_pvt_removedAll(i);
            
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                this->_pvt_added(i, *it, map);
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
bool G2FF::group_setContents(quint32 i, const MoleculeGroup &new_group, 
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
        
        //assert that none of these molecules overlap with the
        //molecules in the other group
        this->assertNoOverlap(molgroup[i!=1], new_group.molecules());
    }
    
    boost::shared_ptr<FF> old_state( this->clone() );
    
    try
    {
        bool changed = molgroup[i].setContents(new_group);
        
        if (not changed)
        {
            for (MoleculeGroup::const_iterator it = new_group.constBegin();
                 it != new_group.constEnd();
                 ++it)
            {
                if (this->_pvt_wouldChangeProperties(i, it->number(), map))
                {
                    changed = true;
                    break;
                }
            }
        }
        
        if (changed)
        {
            this->_pvt_removedAll(i);
            
            for (MoleculeGroup::const_iterator it = new_group.constBegin();
                 it != new_group.constEnd();
                 ++it)
            {
                this->_pvt_added(i, *it, map);
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

/** Return whether or not this forcefield is using any temporary workspace
    that needs to be accepted */
bool G2FF::needsAccepting() const
{
    return molgroup[0].needsAccepting() or molgroup[1].needsAccepting();
}

/** Tell the forcefield that the last move was accepted. This tells the
    forcefield to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void G2FF::accept()
{
    if (molgroup[0].needsAccepting())
        molgroup[0].accept();
    
    if (molgroup[1].needsAccepting())
        molgroup[1].accept();
}
