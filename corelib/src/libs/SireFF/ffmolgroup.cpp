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

#include "ffmolgroup.h"

#include "SireBase/propertymap.h"

#include "SireMol/mgnum.h"
#include "SireMol/mgname.h"
#include "SireMol/molnum.h"
#include "SireMol/molecule.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/viewsofmol.h"
#include "SireMol/molecules.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireMol;
using namespace SireStream;

/////////
///////// Implementation of FFMolGroupPvt
/////////

static const RegisterMetaType<FFMolGroupPvt> r_ffmolgrouppvt( MAGIC_ONLY,
                                                 FFMolGroupPvt::typeName() );
                                                 
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const FFMolGroupPvt &ffmolgrouppvt)
{
    writeHeader(ds, r_ffmolgrouppvt, 1);

    ds << ffmolgrouppvt.mgidx << static_cast<const MoleculeGroup&>(ffmolgrouppvt);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, FFMolGroupPvt &ffmolgrouppvt)
{
    VersionID v = readHeader(ds, r_ffmolgrouppvt);
    
    if (v == 1)
    {
        ds >> ffmolgrouppvt.mgidx >> static_cast<MoleculeGroup&>(ffmolgrouppvt);
        
        //we have to rely on our parent to set the ff pointer...
        ffmolgrouppvt.ffield = 0;
    }
    else
        throw version_error(v, "1", r_ffmolgrouppvt, CODELOC);
        
    return ds;
}

/** Null constructor */
FFMolGroupPvt::FFMolGroupPvt() : MoleculeGroup(), ffield(0)
{}

/** Construct a named group that is a child of the forcefield 'ff' */
FFMolGroupPvt::FFMolGroupPvt(const QString &name, quint32 i, FF *ff)
              : MoleculeGroup(name), mgidx(i), ffield(ff)
{}

/** Copy constructor (our parent forcefield must update our pointer) */
FFMolGroupPvt::FFMolGroupPvt(const FFMolGroupPvt &other)
              : MoleculeGroup(other), mgidx(other.mgidx), ffield(0)
{}

/** Destructor */
FFMolGroupPvt::~FFMolGroupPvt()
{}

/** Copy assignment operator */
FFMolGroupPvt& FFMolGroupPvt::operator=(const FFMolGroupPvt &other)
{
    MoleculeGroup::operator=(other);
    
    //don't change the parent or index - this is because our parent
    //sets this pointer at construction
    
    return *this;
}

/** Set the parent forcefield of this group (this is the forcefield
    that contains this group) */
void FFMolGroupPvt::setParent(FF *new_parent)
{
    ffield = new_parent;
}

/** Set the index of this group in the parent forcefield */
void FFMolGroupPvt::setIndex(quint32 new_idx)
{
    mgidx = MGIdx(new_idx);
}

/** Assert that this is not null (which would be a program bug!)

    \throw SireError::program_bug
*/
void FFMolGroupPvt::assertNotNull() const
{
    if (ffield == 0)
        throw SireError::program_bug( QObject::tr(
            "There is a problem - the parent forcefield of %1 (%2) "
            "is null!").arg(this->name()).arg(this->number()),
                CODELOC );
}

MoleculeGroup* FFMolGroupPvt::clone() const
{
    //return a FFMolGroup, not an FFMolGroupPvt - this
    //allows the FFMolGroup to contain a *copy* of the
    //forcefield
    return new FFMolGroup(*this);
}

/////////
///////// Implementation of FFMolGroup
/////////

static const RegisterMetaType<FFMolGroup> r_ffmolgroup;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const FFMolGroup &ffmolgroup)
{
    writeHeader(ds, r_ffmolgroup, 1);
    
    SharedDataStream sds(ds);
    
    sds << ffmolgroup.mgidx << ffmolgroup.ffield 
        << static_cast<const MoleculeGroup&>(ffmolgroup);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, FFMolGroup &ffmolgroup)
{
    VersionID v = readHeader(ds, r_ffmolgroup);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ffmolgroup.mgidx >> ffmolgroup.ffield 
            >> static_cast<MoleculeGroup&>(ffmolgroup);
    }
    else
        throw version_error(v, "1", r_ffmolgroup, CODELOC);
        
    return ds;
}

/** Null constructor - this creates a useless FFMolGroup. You can
    only construct a valid FFMolGroup by getting a reference to 
    a valid FFMolGroupPvt from a forcefield */
FFMolGroup::FFMolGroup() : ConcreteProperty<FFMolGroup,MoleculeGroup>()
{}

/** Construct from an FFMolGroupPvt - this grabs a copy of the
    forcefield that contains the FFMolGroupPvt */
FFMolGroup::FFMolGroup(const FFMolGroupPvt &ffmolgroup)
           : ConcreteProperty<FFMolGroup,MoleculeGroup>(ffmolgroup), 
             mgidx(ffmolgroup.index()), ffield( ffmolgroup.forceField() )
{}

/** Construct from a MoleculeGroup - you can only construct from a FFMolGroup
    or from an FFMolGroupPvt 
    
    \throw SireError::invalid_arg
*/
FFMolGroup::FFMolGroup(const MoleculeGroup &other)
           : ConcreteProperty<FFMolGroup,MoleculeGroup>(other)
{
    this->operator=(other);
}

/** Copy constructor */
FFMolGroup::FFMolGroup(const FFMolGroup &other)
           : ConcreteProperty<FFMolGroup,MoleculeGroup>(other), 
             mgidx(other.mgidx), ffield(other.ffield)
{}

/** Destructor */
FFMolGroup::~FFMolGroup()
{}

/** Copy assignment operator */
FFMolGroup& FFMolGroup::operator=(const FFMolGroup &other)
{
    MoleculeGroup::operator=(other);
    ffield = other.ffield;
    mgidx = other.mgidx;
    
    return *this;
}

/** Copy assignment from another MoleculeGroup - you can only assign from
    another FFMolGroup, or from an SireFF::detail::FFMolGroupPvt 
    
    \throw SireError::invalid_arg
*/
FFMolGroup& FFMolGroup::operator=(const MoleculeGroup &other)
{
    if (other.isA<FFMolGroupPvt>())
    {
        const FFMolGroupPvt &ff_other = other.asA<FFMolGroupPvt>();
        MoleculeGroup::operator=(other);
        ffield = ff_other.forceField();
        mgidx = ff_other.index();
    }
    else if (other.isA<FFMolGroup>())
    {
        return this->operator=(other.asA<FFMolGroup>());
    }
    else
        throw SireError::invalid_arg( QObject::tr(
            "You can only assign an FFMolGroup from another FFMolGroup "
            "or from an FFMolGroupPvt. You cannot assign from a %1.")
                .arg(other.what()), CODELOC );

    return *this;
}

/** Assert that this is not a null FFMolGroup

    \throw SireError::nullptr_error
*/
void FFMolGroup::assertNotNull() const
{
    if (ffield.constData() == 0)
        throw SireError::nullptr_error( QObject::tr(
            "You cannot use a null FFMolGroup. Please construct a valid "
            "FFMolGroup by copying one from a valid forcefield."), CODELOC );
}

/** Return the forcefield that contains this molecule group */
const FF& FFMolGroup::forceField() const
{
    assertNotNull();
    return ffield.read();
}

/** Return the index of this group in the parent forcefield */
MGIdx FFMolGroup::index() const
{
    return mgidx;
}

/** Update this group so that it is current */
void FFMolGroup::updateGroup()
{
    MoleculeGroup::operator=( ffield->group(this->number()) );
}

/** Set the name of this molecule group */
void FFMolGroup::setName(const QString &name)
{
    assertNotNull();
    ffield.edit().group_setName(mgidx, name);
    updateGroup();
}

void FFMolGroup::add(const MoleculeView &molview, const PropertyMap &map)
{
    assertNotNull();
    
    ffield.edit().add(molview, this->number(), map);
    
    //get the new version of the group
    updateGroup();
}

void FFMolGroup::add(const ViewsOfMol &molviews, const PropertyMap &map)
{
    assertNotNull();
    
    ffield.edit().add(molviews, this->number(), map);
    
    //get the new version of the group
    updateGroup();
}

void FFMolGroup::add(const Molecules &molecules, const PropertyMap &map)
{
    assertNotNull();
    
    ffield.edit().add(molecules, this->number(), map);
    
    //get the new version of the group
    updateGroup();
}

void FFMolGroup::add(const MoleculeGroup &molgroup, const PropertyMap &map)
{
    this->add(molgroup.molecules(), map);
}

void FFMolGroup::add(const MoleculeView &molview)
{
    this->add(molview, PropertyMap());
}

void FFMolGroup::add(const ViewsOfMol &molviews)
{
    this->add(molviews, PropertyMap());
}

void FFMolGroup::add(const Molecules &molecules)
{
    this->add(molecules, PropertyMap());
}

void FFMolGroup::add(const MoleculeGroup &molgroup)
{
    this->add(molgroup, PropertyMap());
}

bool FFMolGroup::addIfUnique(const MoleculeView &molview, const PropertyMap &map)
{
    assertNotNull();
    
    //update to the existing version of the molecule
    PartialMolecule view(molview);
    view.update( ffield->matchToExistingVersion(view.data()) );
    
    //add the molecule
    bool ret = ffield.edit().group_addIfUnique(mgidx, view, map);
    
    //update the index
    ffield.edit().addToIndex( this->number(), view.number() );
    
    //get the new version of the group
    updateGroup();
    
    return ret;
}

ViewsOfMol FFMolGroup::addIfUnique(const ViewsOfMol &molviews, const PropertyMap &map)
{
    assertNotNull();

    //get the current version of the molecule
    ViewsOfMol views(molviews);
    views.update( ffield->matchToExistingVersion(views.data()) );

    //add the views
    ViewsOfMol ret = ffield.edit().group_addIfUnique(mgidx, views, map);

    //update the index
    ffield.edit().addToIndex( this->number(), views.number() );

    //get the latest version of the group
    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::addIfUnique(const Molecules &molecules, 
                                          const PropertyMap &map)
{
    assertNotNull();
    
    //match to the versions already in this forcefield
    Molecules mols = ffield->matchToExistingVersion(molecules);
    
    //add the molecules
    QList<ViewsOfMol> ret = ffield.edit().group_addIfUnique(mgidx, mols, map); 

    //update the index
    ffield.edit().addToIndex(this->number(), mols.molNums());

    //get the latest version of the group
    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::addIfUnique(const MoleculeGroup &molgroup,
                                          const PropertyMap &map)
{
    return this->addIfUnique(molgroup.molecules(), map);
}

bool FFMolGroup::addIfUnique(const MoleculeView &molview)
{
    return this->addIfUnique(molview, PropertyMap());
}

ViewsOfMol FFMolGroup::addIfUnique(const ViewsOfMol &molviews)
{
    return this->addIfUnique(molviews, PropertyMap());
}

QList<ViewsOfMol> FFMolGroup::addIfUnique(const Molecules &molecules)
{
    return this->addIfUnique(molecules, PropertyMap());
}

QList<ViewsOfMol> FFMolGroup::addIfUnique(const MoleculeGroup &molgroup)
{
    return this->addIfUnique(molgroup, PropertyMap());
}

bool FFMolGroup::remove(const MoleculeView &molview)
{
    assertNotNull();

    //remove the molecule
    bool ret = ffield.edit().group_remove(mgidx, molview);

    if (ret)
    {
        if (not ffield->group(this->number()).contains(molview.data().number()))
            //update the index
            ffield.edit().removeFromIndex(this->number(), molview.data().number());
    }

    //get the latest version of the group
    updateGroup();
    
    return ret;
}

ViewsOfMol FFMolGroup::remove(const ViewsOfMol &molviews)
{
    assertNotNull();
    
    //remove the views
    ViewsOfMol ret = ffield.edit().group_remove(mgidx, molviews);
    
    if (not ret.isEmpty())
    {
        if (not ffield->group(this->number()).contains(molviews.number()))
            //update the index
            ffield.edit().removeFromIndex(this->number(), molviews.number());
    }
    
    //get the latest version of the group
    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::remove(const Molecules &molecules)
{
    assertNotNull();
    
    //remove the molecules
    QList<ViewsOfMol> ret = ffield.edit().group_remove(mgidx, molecules);
    
    if (not ret.isEmpty())
    {
        const MoleculeGroup &group = ffield->group(this->number());
        QSet<MolNum> removed_molnums;
        
        foreach (const ViewsOfMol &removed_mol, ret)
        {
            if (not group.contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }
        
        if (not removed_molnums.isEmpty())
            ffield.edit().removeFromIndex(this->number(), removed_molnums);
    }
    
    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::remove(const MoleculeGroup &molgroup)
{
    return this->remove(molgroup.molecules());
}

bool FFMolGroup::removeAll(const MoleculeView &molview)
{
    assertNotNull();

    //remove the molecule
    bool ret = ffield.edit().group_removeAll(mgidx, molview);

    if (ret)
    {
        if (not ffield->group(this->number()).contains(molview.data().number()))
            //update the index
            ffield.edit().removeFromIndex(this->number(), molview.data().number());
    }

    //get the latest version of the group
    updateGroup();
    
    return ret;
}

ViewsOfMol FFMolGroup::removeAll(const ViewsOfMol &molviews)
{
    assertNotNull();
    
    //remove the views
    ViewsOfMol ret = ffield.edit().group_removeAll(mgidx, molviews);
    
    if (not ret.isEmpty())
    {
        if (not ffield->group(this->number()).contains(molviews.number()))
            //update the index
            ffield.edit().removeFromIndex(this->number(), molviews.number());
    }
    
    //get the latest version of the group
    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::removeAll(const Molecules &molecules)
{
    assertNotNull();
    
    //remove the molecules
    QList<ViewsOfMol> ret = ffield.edit().group_removeAll(mgidx, molecules);
    
    if (not ret.isEmpty())
    {
        const MoleculeGroup &group = ffield->group(this->number());
        QSet<MolNum> removed_molnums;
        
        foreach (const ViewsOfMol &removed_mol, ret)
        {
            if (not group.contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }
        
        if (not removed_molnums.isEmpty())
            ffield.edit().removeFromIndex(this->number(), removed_molnums);
    }
    
    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::removeAll(const MoleculeGroup &molgroup)
{
    return this->removeAll(molgroup.molecules());
}

ViewsOfMol FFMolGroup::remove(MolNum molnum)
{
    assertNotNull();

    ViewsOfMol ret = ffield.edit().group_remove(mgidx, molnum);

    if (not ret.isEmpty())
    {
        ffield.edit().removeFromIndex(this->number(), molnum);
    }

    updateGroup();
    
    return ret;
}

QList<ViewsOfMol> FFMolGroup::remove(const QSet<MolNum> &molnums)
{
    assertNotNull();

    QList<ViewsOfMol> ret = ffield.edit().group_remove(mgidx, molnums);

    QSet<MolNum> removed_molnums;
    
    foreach (const ViewsOfMol &removed_mol, ret)
    {
        removed_molnums.insert(removed_mol.number());
    }
    
    if (not removed_molnums.isEmpty())
        ffield.edit().removeFromIndex(this->number(), removed_molnums);

    updateGroup();
    
    return ret;
}

void FFMolGroup::removeAll()
{
    assertNotNull();
    ffield.edit().removeAll(this->number());
    updateGroup();
}

bool FFMolGroup::update(const MoleculeData &moldata, bool auto_commit)
{
    assertNotNull();
    
    //update in a copy of the group (so get the return value)
    ForceField copy = ffield;
    bool ret = copy.edit().group_update(mgidx, moldata, auto_commit);

    if (ret)
    {
        //we have to do the update for real
        ffield.edit().update(moldata, auto_commit);
        updateGroup();
    }
    
    return false;
}

QList<Molecule> FFMolGroup::update(const Molecules &molecules, bool auto_commit)
{
    assertNotNull();

    //perform the update in a copy so that we can get the return 
    //value
    ForceField copy = ffield;
    QList<Molecule> ret = copy.edit().group_update(mgidx, molecules, auto_commit);
    
    if (not ret.isEmpty())
    {
        //we have to update the entire forcefield
        ffield.edit().update(molecules, auto_commit);
        updateGroup();
    }
    
    return ret;
}

QList<Molecule> FFMolGroup::update(const MoleculeGroup &molgroup, bool auto_commit)
{
    assertNotNull();

    //perform the update in a copy so that we can get the return 
    //value
    ForceField copy = ffield;
    QList<Molecule> ret = copy.edit().group_update(mgidx, molgroup, auto_commit);
    
    if (not ret.isEmpty())
    {
        //we have to update the entire forcefield
        ffield.edit().update(molgroup, auto_commit);
        updateGroup();
    }
    
    return ret;
}

bool FFMolGroup::setContents(const MoleculeView &molview, const PropertyMap &map)
{
    assertNotNull();

    //what's the old version of this group?
    quint32 majver = this->majorVersion();
    quint32 minver = this->minorVersion();
    
    //make the change
    ffield.edit().setContents(this->number(), molview, map);
    
    //update this group
    updateGroup();
    
    //has there been a change
    return majver != this->majorVersion() or 
           minver != this->minorVersion();
}

bool FFMolGroup::setContents(const ViewsOfMol &molviews, const PropertyMap &map)
{
    assertNotNull();

    //what's the old version of this group?
    quint32 majver = this->majorVersion();
    quint32 minver = this->minorVersion();
    
    //make the change
    ffield.edit().setContents(this->number(), molviews, map);
    
    //update this group
    updateGroup();
    
    //has there been a change
    return majver != this->majorVersion() or 
           minver != this->minorVersion();
}

bool FFMolGroup::setContents(const Molecules &molecules, const PropertyMap &map)
{
    assertNotNull();

    //what's the old version of this group?
    quint32 majver = this->majorVersion();
    quint32 minver = this->minorVersion();
    
    //make the change
    ffield.edit().setContents(this->number(), molecules, map);
    
    //update this group
    updateGroup();
    
    //has there been a change
    return majver != this->majorVersion() or 
           minver != this->minorVersion();
}

bool FFMolGroup::setContents(const MoleculeGroup &molgroup, const PropertyMap &map)
{
    assertNotNull();

    //what's the old version of this group?
    quint32 majver = this->majorVersion();
    quint32 minver = this->minorVersion();
    
    //make the change
    ffield.edit().setContents(this->number(), molgroup, map);
    
    //update this group
    updateGroup();
    
    //has there been a change
    return majver != this->majorVersion() or 
           minver != this->minorVersion();
}

bool FFMolGroup::setContents(const MoleculeView &molview)
{
    return this->setContents(molview, PropertyMap());
}

bool FFMolGroup::setContents(const ViewsOfMol &molviews)
{
    return this->setContents(molviews, PropertyMap());
}

bool FFMolGroup::setContents(const Molecules &molecules)
{
    return this->setContents(molecules, PropertyMap());
}

bool FFMolGroup::setContents(const MoleculeGroup &molgroup)
{
    return this->setContents(molgroup, PropertyMap());
}

const char* FFMolGroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FFMolGroup>() );
}

bool FFMolGroup::needsAccepting() const
{
    return MoleculeGroup::needsAccepting() or ffield.read().needsAccepting();
}

void FFMolGroup::accept()
{
    MoleculeGroup::accept();
    
    if (ffield.read().needsAccepting())
        ffield.edit().accept();
}

FFMolGroup* FFMolGroup::clone() const
{
    return new FFMolGroup(*this);
}

