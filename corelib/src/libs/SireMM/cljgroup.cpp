/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "cljgroup.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QElapsedTimer>

using namespace SireMM;
using namespace SireMol;
using namespace SireStream;
using namespace SireBase;

static const RegisterMetaType<CLJGroup> r_group(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJGroup &group)
{
    writeHeader(ds, r_group, 2);
    
    SharedDataStream sds(ds);
    
    sds << group.cljexts << group.cljboxes
        << group.cljworkspace << group.changed_mols
        << group.props << qint32(group.id_source) << qint32(group.extract_source);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJGroup &group)
{
    VersionID v = readHeader(ds, r_group);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        qint32 id_source;
        qint32 extract_source;
        
        sds >> group.cljexts >> group.cljboxes
            >> group.cljworkspace >> group.changed_mols
            >> group.props >> id_source >> extract_source;
        
        group.id_source = CLJAtoms::ID_SOURCE(id_source);
        group.extract_source = CLJExtractor::EXTRACT_SOURCE(extract_source);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        qint32 id_source;
        bool split_by_residue;
        
        sds >> group.cljexts >> group.cljboxes
            >> group.cljworkspace >> group.changed_mols
            >> group.props >> id_source >> split_by_residue;
        
        group.id_source = CLJAtoms::ID_SOURCE(id_source);
        
        if (split_by_residue)
        {
            group.extract_source = CLJExtractor::EXTRACT_BY_RESIDUE;
        }
        else
        {
            group.extract_source = CLJExtractor::EXTRACT_BY_MOLECULE;
        }
    }
    else
        throw version_error(v, "1", r_group, CODELOC);
    
    return ds;
}

/** Null constructor */
CLJGroup::CLJGroup()
         : id_source(CLJAtoms::USE_MOLNUM), extract_source(CLJExtractor::EXTRACT_BY_CUTGROUP)
{
    cljworkspace.mustRecalculateFromScratch(cljboxes);
}

/** Construct, suppling the name of the molecule group and the source of the
    CLJAtoms ID_SOURCE property (e.g. USE_MOLNUM for intermolecular forcefields or
    USE_ATOMNUM for intramolecular forcefields) */
CLJGroup::CLJGroup(CLJAtoms::ID_SOURCE ids)
         : id_source(ids), extract_source(CLJExtractor::EXTRACT_BY_CUTGROUP)
{
    cljworkspace.mustRecalculateFromScratch(cljboxes);
}

/** Construct specifying how the atoms will be extracted from the molecule */
CLJGroup::CLJGroup(CLJExtractor::EXTRACT_SOURCE ext)
         : id_source(CLJAtoms::USE_MOLNUM), extract_source(ext)
{
    cljworkspace.mustRecalculateFromScratch(cljboxes);
}

/** Construct, suppling the name of the molecule group and the source of the
    CLJAtoms ID_SOURCE property (e.g. USE_MOLNUM for intermolecular forcefields or
    USE_ATOMNUM for intramolecular forcefields), and also specifying if we are going
    to split molecules by cutgroup, residue or molecule, (normally best to
    extract by cutgroup unless you know that all molecules are going to be small) */
CLJGroup::CLJGroup(CLJAtoms::ID_SOURCE ids, CLJExtractor::EXTRACT_SOURCE ext)
         : id_source(ids), extract_source(ext)
{
    cljworkspace.mustRecalculateFromScratch(cljboxes);
}

/** Copy constructor */
CLJGroup::CLJGroup(const CLJGroup &other)
         : cljexts(other.cljexts), cljboxes(other.cljboxes),
           cljworkspace(other.cljworkspace), changed_mols(other.changed_mols),
           props(other.props), id_source(other.id_source),
           extract_source(other.extract_source)
{}

/** Destructor */
CLJGroup::~CLJGroup()
{}

/** Copy assignment operator */
CLJGroup& CLJGroup::operator=(const CLJGroup &other)
{
    if (this != &other)
    {
        cljexts = other.cljexts;
        cljboxes = other.cljboxes;
        cljworkspace = other.cljworkspace;
        changed_mols = other.changed_mols;
        props = other.props;
        id_source = other.id_source;
        extract_source = other.extract_source;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJGroup::operator==(const CLJGroup &other) const
{
    return (this == &other) or
           (cljexts == other.cljexts and
            cljboxes == other.cljboxes and cljworkspace == other.cljworkspace and
            changed_mols == other.changed_mols and props == other.props and
            id_source == other.id_source and extract_source == other.extract_source);
}

/** Comparison operator */
bool CLJGroup::operator!=(const CLJGroup &other) const
{
    return not operator==(other);
}

const char* CLJGroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJGroup>() );
}

const char* CLJGroup::what() const
{
    return CLJGroup::typeName();
}

QString CLJGroup::toString() const
{
    return QObject::tr("CLJGroup( needsAccepting() == %1 )")
                .arg(needsAccepting());
}

/** Return the current version of all of the molecules in this group */
Molecules CLJGroup::molecules() const
{
    Molecules mols;

    for (QHash<MolNum,CLJExtractor>::const_iterator it = changed_mols.constBegin();
         it != changed_mols.constEnd(); ++it)
    {
        mols.add( it.value().newMolecule() );
    }
    
    for (ChunkedHash<MolNum,CLJExtractor>::const_iterator it = cljexts.constBegin();
         it != cljexts.constEnd(); ++it)
    {
        if (not mols.contains(it.key()))
        {
            mols.add( it.value().newMolecule() );
        }
    }
    
    return mols;
}

/** Return the size of the box used by CLJBoxes to partition space */
Length CLJGroup::boxLength() const
{
    return cljboxes.length();
}

/** Set the size of the box used by CLJBoxes to partition space */
void CLJGroup::setBoxLength(Length box_length)
{
    if (box_length != boxLength())
    {
        cljboxes = CLJBoxes(box_length);
        this->mustReallyRecalculateFromScratch();
    }
}

/** Return the property map used for the molecule with number 'molnum' */
PropertyMap CLJGroup::mapForMolecule(MolNum molnum) const
{
    return props.value(molnum, PropertyMap());
}

/** Add the passed molecule to this group */
void CLJGroup::add(const MoleculeView &molview, const PropertyMap &map)
{
    MolNum molnum = molview.data().number();
    
    bool must_add_map = false;
    
    if (cljexts.contains(molnum) or changed_mols.contains(molnum))
    {
        //we already contain this molecule - we should check that there is no
        //change in property map
        if (props.value(molnum,PropertyMap()) != map)
        {
            throw SireError::incompatible_error( QObject::tr(
                    "Cannot add the molecule %1 because the passed PropertyMap (%2) "
                    "is not the same as that used the last time the molecule was "
                    "added (%3)")
                        .arg(molview.toString())
                        .arg(map.toString())
                        .arg(props.value(molnum,PropertyMap()).toString()), CODELOC );
        }
    }
    else if (map != PropertyMap())
    {
        must_add_map = true;
    }

    if (changed_mols.contains(molnum))
    {
        changed_mols[molnum].add(molview, cljboxes, cljworkspace);
    }
    else if (cljexts.contains(molnum))
    {
        if (cljworkspace.recalculatingFromScratch())
        {
            cljexts[molnum].add(molview, cljboxes, cljworkspace);
        }
        else
        {
            CLJExtractor changed = cljexts.value(molnum);
            changed.add(molview, cljboxes, cljworkspace);
            changed_mols.insert(molnum, changed);
        }
    }
    else
    {
        CLJExtractor changed(molview, id_source, extract_source, map);
        
        if (cljworkspace.recalculatingFromScratch())
        {
            changed.commit(cljboxes, cljworkspace);
            cljexts.insert(molnum, changed);
        }
        else
        {
            changed_mols.insert(molnum, changed);
        }
        
        if (must_add_map)
            props.insert(molnum, map);
    }
}

/** Add all of the passed molecules to this group */
void CLJGroup::add(const Molecules &molecules, const PropertyMap &map)
{
    CLJGroup old_group(*this);
    
    try
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            this->add(it.value(), map);
        }
    }
    catch(...)
    {
        this->operator=(old_group);
        throw;
    }
}

/** Add all of the passed molecules to this group */
void CLJGroup::add(const MoleculeGroup &molecules, const PropertyMap &map)
{
    this->add(molecules.molecules(), map);
}

/** Update the molecule in this group to match the version of the molecule 'molview' */
void CLJGroup::update(const MoleculeView &molview)
{
    MolNum molnum = molview.data().number();
    
    if (cljexts.contains(molnum) or changed_mols.contains(molnum))
    {
        if (cljexts.count() == 1)
        {
            const CLJExtractor &ext = cljexts.constBegin().value();
            
            if (ext.extractingByMolecule() or
                (ext.extractingByResidue() and molview.data().info().nResidues() == 1) or
                (ext.extractingByCutGroup() and molview.data().info().nCutGroups() == 1))
            {
                //there is only a single CLJAtoms group in this object, so we may
                //as well rebuild it all from scratch
                this->mustRecalculateFromScratch();
                cljexts[molnum].update(molview, cljboxes, cljworkspace);
                return;
            }
        }

        //we must update this molecule
        if (changed_mols.contains(molnum))
        {
            changed_mols[molnum].update(molview, cljboxes, cljworkspace);
        }
        else if (cljexts.contains(molnum))
        {
            if (cljworkspace.recalculatingFromScratch())
            {
                cljexts[molnum].update(molview, cljboxes, cljworkspace);
            }
            else
            {
                CLJExtractor changed = cljexts.value(molnum);
                changed.update(molview, cljboxes, cljworkspace);
                
                if (changed.hasChangedAtoms())
                {
                    changed_mols.insert(molnum,changed);
                }
                else
                {
                    //this change hasn't actually affected the atoms
                    changed.commit(cljboxes, cljworkspace);
                    cljexts[molnum] = changed;
                }
            }
        }
    }
}

/** Update all of the molecules in this group that are in 'molecules' to match
    the version held in 'molecules' */
void CLJGroup::update(const Molecules &molecules)
{
    CLJGroup old_group(*this);
    
    try
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd(); ++it)
        {
            this->update(it.value());
        }
    }
    catch(...)
    {
        this->operator=(old_group);
        throw;
    }
}

/** Update all of the molecules in this group that are in 'molecules' to match
    the version held in 'molecules' */
void CLJGroup::update(const MoleculeGroup &molecules)
{
    this->update(molecules.molecules());
}

/** Remove the molecule view 'molview' */
void CLJGroup::remove(const MoleculeView &molview)
{
    MolNum molnum = molview.data().number();
    
    PartialMolecule newmol;
    
    if (changed_mols.contains(molnum))
    {
        changed_mols[molnum].remove(molview, cljboxes, cljworkspace);
        newmol = changed_mols[molnum].newMolecule();
    }
    else if (cljexts.contains(molnum))
    {
        if (cljworkspace.recalculatingFromScratch())
        {
            cljexts[molnum].remove(molview, cljboxes, cljworkspace);
            newmol = cljexts[molnum].newMolecule();
            
            if (newmol.selection().nSelected() == 0)
                cljexts.remove(molnum);
        }
        else
        {
            CLJExtractor changed = cljexts.value(molnum);
            changed.remove(molview, cljboxes, cljworkspace);
            changed_mols.insert(molnum, changed);
            newmol = changed.newMolecule();
        }
    }

    if (newmol.nAtoms() == 0)
    {
        props.remove(molnum);
    }
}

/** Remove all of the molecules in 'molecules' from this group */
void CLJGroup::remove(const Molecules &molecules)
{
    CLJGroup old_group(*this);
    
    try
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd(); ++it)
        {
            this->remove(it.value());
        }
    }
    catch(...)
    {
        this->operator=(old_group);
        throw;
    }
}

/** Remove all of the molecules in 'molecules' from this group */
void CLJGroup::remove(const MoleculeGroup &molecules)
{
    this->remove(molecules.molecules());
}

/** Remove the molecule with number 'molnum' from this group */
void CLJGroup::remove(MolNum molnum)
{
    if (changed_mols.contains(molnum))
    {
        changed_mols[molnum].removeAll(cljboxes, cljworkspace);
    }
    else if (cljexts.contains(molnum))
    {
        if (cljworkspace.recalculatingFromScratch())
        {
            cljexts[molnum].removeAll(cljboxes, cljworkspace);
            cljexts.remove(molnum);
        }
        else
        {
            CLJExtractor changed = cljexts.value(molnum);
            changed.removeAll(cljboxes, cljworkspace);
            changed_mols.insert(molnum, changed);
        }
    }
    
    if (props.contains(molnum))
    {
        props.remove(molnum);
    }
}

/** Remove all molecules from this group */
void CLJGroup::removeAll()
{
    cljexts.clear();
    changed_mols.clear();
    cljboxes = CLJBoxes(cljboxes.length());
    cljworkspace = CLJWorkspace();
    cljworkspace.mustRecalculateFromScratch(cljboxes);
    props.clear();
}

/** Return the number of molecules that have changed since the 
    last time this CLJGroup was accepted */
int CLJGroup::nChangedMolecules() const
{
    if (recalculatingFromScratch())
    {
        //assume that everything has changed
        return cljexts.count();
    }
    else
    {
        return changed_mols.count();
    }
}

/** Return all of the molecules that have changed since the last time
    this CLJGroup was accepted */
Molecules CLJGroup::changedMolecules() const
{
    if (recalculatingFromScratch())
    {
        return this->molecules();
    }
    else
    {
        Molecules mols;
        
        for (QHash<MolNum,CLJExtractor>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd(); ++it)
        {
            mols.add( it.value().newMolecule() );
        }
        
        return mols;
    }
}

/** Return whether or not this group needs to be accepted */
bool CLJGroup::needsAccepting() const
{
    if (cljworkspace.needsAccepting() or (not changed_mols.isEmpty()))
        return true;
    
    else
    {
        for (ChunkedHash<MolNum,CLJExtractor>::const_iterator it = cljexts.constBegin();
             it != cljexts.constEnd(); ++it)
        {
            if (it.value().needsCommitting())
                return true;
        }
        
        return false;
    }
}

/** Accept all of the changes in the group. This will ensure that 
    all deltas have been removed and all of the atoms are correctly
    added to the CLJBoxes boxes */
void CLJGroup::accept()
{
    QSet<MolNum> removed_mols;

    if (cljworkspace.recalculatingFromScratch())
    {
        //this adds any molecules to the box that have not yet been
        //added to the box
        for (ChunkedHash<MolNum,CLJExtractor>::iterator it = cljexts.begin();
             it != cljexts.end();
             ++it)
        {
            it.value().commit(cljboxes, cljworkspace);
            
            if (it.value().isEmpty())
                removed_mols.insert(it.key());
        }

        cljworkspace.accept(cljboxes);
        changed_mols.clear();
    }
    else if (not changed_mols.isEmpty())
    {
        for (QHash<MolNum,CLJExtractor>::iterator it = changed_mols.begin();
             it != changed_mols.end(); ++it)
        {
            it.value().commit(cljboxes, cljworkspace);
            cljexts[it.key()] = it.value();
            
            if (it.value().isEmpty())
                removed_mols.insert(it.key());
        }

        cljworkspace.accept(cljboxes);
        changed_mols.clear();
    }
    else if (cljworkspace.needsAccepting())
    {
        cljworkspace.accept(cljboxes);
    }

    if (not removed_mols.isEmpty())
    {
        foreach (MolNum removed_mol, removed_mols)
        {
            cljexts.remove(removed_mol);
            props.remove(removed_mol);
        }
        
        if (cljexts.isEmpty())
        {
            cljexts.clear();
            changed_mols.clear();
            props.clear();
        }
    }
}

/** Tell this CLJGroup that atoms in a connected group have been
    updated. This tells the CLJGroup to remove some of the caching
    that is used to improve performance, as this caching could cause
    energy errors when using deltas generated from other CLJGroups */
void CLJGroup::updatedConnectedGroup()
{
    cljworkspace.removeSameIDAtoms(cljboxes);
}

/** Return whether or not the changes since the last time 'accept()'
    was called all change parts of a molecule that all have the same
    ID. If they do, then we only need to use 'changedAtoms()' in the
    delta energy calculation as the changed atoms don't interact with
    each other. Otherwise, we need 'mergeChanges()' to get the changedAtoms
    together with oldAtoms and newAtoms, so that we can calculate the
    change in energy within changedAtoms() itself as well */
bool CLJGroup::isSingleIDChange() const
{
    if (cljworkspace.recalculatingFromScratch())
        return false;
    else
        return cljworkspace.isSingleID();
}

/** Return a tuple of (changedAtoms(),oldAtoms(),newAtoms()). This is
    needed if more than a single ID group has changed and thus we
    need to calculate the change in interaction within changedAtoms */
tuple<CLJAtoms,CLJAtoms,CLJAtoms> CLJGroup::mergeChanges() const
{
    if (cljworkspace.recalculatingFromScratch())
        return tuple<CLJAtoms,CLJAtoms,CLJAtoms>();
    else
        return cljworkspace.merge();
}

/** Return the set of all atoms that have changed since the last 
    time 'accept()' was called. This will return an empty set if 
    the workspace was told to "mustRecalculateFromScratch()" */
CLJAtoms CLJGroup::changedAtoms() const
{
    if (cljworkspace.recalculatingFromScratch())
        return CLJAtoms();
    else
        return cljworkspace.changedAtoms();
}

/** Return the set of all of the new atoms that have changed
    since the last time 'accept()' was called. This, plus the 
    negative of oldAtoms() will equal changedAtoms() */
CLJAtoms CLJGroup::newAtoms() const
{
    if (cljworkspace.recalculatingFromScratch())
        return CLJAtoms();
    else
        return cljworkspace.newAtoms();
}

/** Return the set of all of the old atoms that have changed
    since the last time 'accept()' was called. The negative
    of this plus newAtoms() will equal changedAtoms() */
CLJAtoms CLJGroup::oldAtoms() const
{
    if (cljworkspace.recalculatingFromScratch())
        return CLJAtoms();
    else
        return cljworkspace.oldAtoms();
}

/** Tell the group that calculations will be made completely from scratch,
    so there is no need to maintain a delta */
void CLJGroup::mustRecalculateFromScratch()
{
    this->accept();
    cljworkspace.mustRecalculateFromScratch(cljboxes);
}

/** Tell the group that calculations will be made completely from scratch,
    and to also re-extract all of the molecules */
void CLJGroup::mustReallyRecalculateFromScratch()
{
    this->mustRecalculateFromScratch();

    cljworkspace = CLJWorkspace();
    cljboxes = CLJBoxes(this->boxLength());

    changed_mols.clear();
    
    Molecules mols;
    
    for (ChunkedHash<MolNum,CLJExtractor>::const_iterator it = cljexts.constBegin();
         it != cljexts.constEnd(); ++it)
    {
        mols.add( it.value().newMolecule() );
    }

    cljexts.clear();

    for (MoleculeGroup::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        CLJExtractor cljmol( it.value(), id_source, extract_source,
                             props.value(it.key(),PropertyMap()) );
    
        cljmol.commit(cljboxes, cljworkspace);;
    
        cljexts.insert( it.key(), cljmol );
    }
}
