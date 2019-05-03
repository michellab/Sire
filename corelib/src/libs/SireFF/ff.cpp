/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "ff.h"
#include "ffcomponent.h"
#include "forcefield.h"

#include "tostring.h"

#include "SireMol/mgnum.h"
#include "SireMol/mgidx.h"
#include "SireMol/mgname.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/viewsofmol.h"
#include "SireMol/molecules.h"
#include "SireMol/moleculegroup.h"

#include "SireMol/mover.hpp"

#include "SireFF/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>
#include <QDebug>

using namespace SireFF;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

void SireFF::detail::throwForceFieldRestoreBug(const char *this_what, 
                                                             const char *ffield_what)
{
    throw SireError::program_bug( QObject::tr(
        "Something went wrong, as we are trying to restore a forcefield "
        "of type %1 using a saved state forcefield of type %2...??")
            .arg(this_what).arg(ffield_what), CODELOC );
}

static const RegisterMetaType<FF> r_ff( MAGIC_ONLY, FF::typeName() );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const FF &ff)
{
    writeHeader(ds, r_ff, 1);
    
    SharedDataStream sds(ds);
    
    sds << ff.uid << ff.versn << ff.ffname << ff.nrg_components
        << ff.isdirty
        << static_cast<const MolGroupsBase&>(ff);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, FF &ff)
{
    VersionID v = readHeader(ds, r_ff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> ff.uid >> ff.versn >> ff.ffname >> ff.nrg_components
            >> ff.isdirty
            >> static_cast<MolGroupsBase&>(ff);
    }
    else
        throw version_error(v, "1", r_ff, CODELOC);
        
    return ds;
}

/** Constructor */
FF::FF() : MolGroupsBase(), versn(0), 
                            version_ptr( new Incremint(0) ),
                            ffname( QObject::tr("unnamed") ),
                            isdirty(true)
{
    uid = QUuid::createUuid();
}

/** Construct a forcefield, and also give it a name */
FF::FF(const QString &name) : MolGroupsBase(), 
                              versn(0), 
                              version_ptr( new Incremint(0) ),
                              ffname(name),
                              isdirty(true)
{
    uid = QUuid::createUuid();
}

/** Copy constructor */
FF::FF(const FF &other) : MolGroupsBase(other),
                          uid(other.uid), versn(other.versn),
                          version_ptr(other.version_ptr),
                          ffname(other.ffname),
                          nrg_components(other.nrg_components),
                          isdirty(other.isdirty)
{}

/** Destructor */
FF::~FF()
{}

/** Copy assignment operator */
FF& FF::operator=(const FF &other)
{
    if (this != &other)
    {
        uid = other.uid;
        versn = other.versn;
        version_ptr = other.version_ptr;
        ffname = other.ffname;
        nrg_components = other.nrg_components;
        isdirty = other.isdirty;
        
        MolGroupsBase::operator=(other);
    }
    
    return *this;
}

/** Comparison operator - two forcefields are only identical if
    they have the same UID and version */
bool FF::operator==(const FF &other) const
{
    return uid == other.uid and versn == other.versn;
}

/** Comparison operator - two forcefields are only identical if
    they have the same UID and version */
bool FF::operator!=(const FF &other) const
{
    return uid != other.uid or versn != other.versn;
}

/** Return the current values of the energy components - note that
    these will be invalid if this->isDirty() is true */
const Values& FF::currentEnergies() const
{
    return nrg_components;
}

/** Return the unique ID for this forcefield */
const QUuid& FF::UID() const
{
    return uid;
}

/** Return the version number of this copy of the forcefield */
quint64 FF::version() const
{
    return versn;
}

/** Return a string representation of this forcefield */
QString FF::toString() const
{
    return QObject::tr("%1(\"%2\", version=%3)")
                    .arg(this->what())
                    .arg(this->name())
                    .arg(this->version());
}

/** Increment the version number of this forcefield */
void FF::incrementVersion()
{
    versn = version_ptr->increment();
}

/** Return the name of this forcefield */
const FFName& FF::name() const
{
    return ffname;
}

/** Set the name of this forcefield */
void FF::setName(const QString &name)
{
    FFName old_name = ffname;

    try
    {
        if (ffname != FFName(name))
        {
            ffname = FFName(name);
            
            uid = QUuid::createUuid();
            version_ptr.reset( new Incremint() );
            
            this->_pvt_updateName();

            this->incrementVersion();
        }
    }
    catch(...)
    {
        ffname = old_name;
        throw;
    }
}

/** Set the name of the forcefield groups identified by 'mgid' 

    \throw SireMol::missing_group
*/
void FF::setName(const MGID &mgid, const QString &name)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    foreach (MGNum mgnum, mgnums)
    {
        MGName old_name = this->group(mgnum).name();
        MGName new_name = MGName(name);
        
        if (old_name != new_name)
        {
            this->changeNameIndex(mgnum, old_name, new_name);
            this->group_setName(this->mgIdx(mgnum), new_name);
        }
    }
}

/** Return the energy of the forcefield component represented
    by the passed symbol
    
    \throw SireFF::missing_component
*/
SireUnits::Dimension::MolarEnergy FF::energy(const Symbol &component)
{
    if (this->isDirty())
    {
        this->recalculateEnergy();
    }
                  
    if (not nrg_components.values().contains(component.ID()))
        throw SireFF::missing_component( QObject::tr(
            "There is no component in this forcefield represented by "
            "the symbol %1. Available components are %2.")
                .arg(component.toString())
                .arg( Sire::toString(this->components().symbols()) ), CODELOC );

    return SireUnits::Dimension::MolarEnergy( nrg_components.value(component) );
}

/** Return the energy of this forcefield in its current state */
SireUnits::Dimension::MolarEnergy FF::energy()
{
    return this->energy( this->components().total() );
}

/** Return the values of the specified energy components */
Values FF::energies(const QSet<Symbol> &components)
{
    if (this->isDirty())
        this->recalculateEnergy();
        
    Values vals;
    
    foreach (const Symbol &component, components)
    {
        if (not nrg_components.values().contains(component.ID()))
            throw SireFF::missing_component( QObject::tr(
                "There is no component in this forcefield represented by "
                "the symbol %1. Available components are %2.")
                    .arg(component.toString())
                    .arg( Sire::toString( nrg_components.values().keys() ) ), CODELOC );
    
        vals.set( component, nrg_components.value(component) );
    }
    
    return vals;
}

/** Return the values of all of the energy components of this forcefield */
Values FF::energies()
{
    if (this->isDirty())
        this->recalculateEnergy();
        
    return nrg_components;
}

/** Add the molecule viewed in 'molview' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters
  
    Note that if this molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version of 'molview'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const MoleculeView &molview, const MGID &mgid,
             const PropertyMap &map)
{
    //get the numbers of the forcefield groups...
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;
    
    //update the molecule so that it is at the same version as 
    //any existing copies
    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );
                
    if (mgnums.count() == 1)
    {
        //no need to checkpoint as the order of operations can
        //be used to preserve state
        MGIdx mgidx = this->mgIdx( *(mgnums.constBegin()) );
        
        this->group_add(mgidx, view, map);
        this->addToIndex( *(mgnums.constBegin()), view.number() );
    }
    else
    {
        //we need to save state as an exception could be thrown
        //when adding to the nth group
        boost::shared_ptr<FF> old_state( this->clone() );
        
        try
        {
            //add the molecule to each group in turn...
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                this->group_add(mgidx, view, map);
                this->addToIndex(mgnum, view.number());
            }
        }
        catch(...)
        {
            //restore the old state of the forcefield
            this->copy(*old_state);
            throw;
        }
    }
}

/** Add the views of the molecule viewed in 'molviews' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters
  
    Note that if this molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version of 'molview'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const ViewsOfMol &molviews, const MGID &mgid,
             const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;
        
    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        this->group_add(mgidx, views, map);
        this->addToIndex(mgnum, views.number());
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                this->group_add(mgidx, views, map);
                this->addToIndex(mgnum, views.number());
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}

/** Add all of the molecules in 'molecules' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters
  
    Note that if a molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version in 'molecules'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const Molecules &molecules, const MGID &mgid,
             const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;
    
    Molecules mols = this->matchToExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();
    
    if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        this->group_add(mgidx, molecules, map);
        this->addToIndex(mgnum, molnums);
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                this->group_add(mgidx, molecules, map);
                this->addToIndex(mgnum, molnums);
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}
             
/** Add all of the molecules in the group 'molgroup' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters
  
    Note that if a molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version in 'molecules'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const MoleculeGroup &molgroup, const MGID &mgid,
             const PropertyMap &map)
{
    this->add(molgroup.molecules(), mgid, map);
}

/** Add the molecule viewed in 'molview' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters

    This only adds the view to groups that don't already contain it.
  
    Note that if this molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version of 'molview'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const MoleculeView &molview, const MGID &mgid,
                     const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;
        
    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );
    
    if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        this->group_addIfUnique(mgidx, view, map);
        this->addToIndex(mgnum, view.number());
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                this->group_addIfUnique(mgidx, view, map);
                this->addToIndex(mgnum, view.number());
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}

/** Add the views of the molecule viewed in 'molviews' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters

    This only adds views to groups that don't already contain them.
  
    Note that if this molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version of 'molview'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid,
                     const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;
        
    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        this->group_addIfUnique(mgidx, views, map);
        this->addToIndex(mgnum, views.number());
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                this->group_addIfUnique(mgidx, views, map);
                this->addToIndex(mgnum, views.number());
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}

/** Add all of the molecules in 'molecules' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters

    This only adds molecules to groups that don't already contain them.

    Note that if a molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version in 'molecules'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const Molecules &molecules, const MGID &mgid,
                     const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;
        
    Molecules mols = this->matchToExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();
    
    if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        this->group_addIfUnique(mgidx, mols, map);
        this->addToIndex(mgnum, molnums);
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                this->group_addIfUnique(mgidx, mols, map);
                this->addToIndex(mgnum, molnums);
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}

/** Add all of the molecules in the group 'molgroup' to the forcefield groups
    identified by 'mgid' using the supplied map to find the properties
    of the molecule that contain the forcefield parameters

    This only adds molecules to groups that don't already contain them.
  
    Note that if a molecule exists already in this forcefield, then
    the version of the molecule that is in this forcefield will be used,
    not the version in 'molecules'
      
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid,
                     const PropertyMap &map)
{
    this->addIfUnique(molgroup.molecules(), mgid, map);
}

/** Completely remove all of the molecules from the groups identified
    by 'mgid'
    
    \throw SireMol::missing_group
*/
bool FF::removeAll(const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    bool mols_removed = false;
    
    foreach (MGNum mgnum, mgnums)
    {
        if (not mols_removed)
        {
            if (not this->at(mgnum).isEmpty())
                mols_removed = true;
        }
        
        this->group_removeAll( this->mgIdx(mgnum) );
        this->clearIndex(mgnum);
    }
    
    return mols_removed;
}

/** Remove the view 'molview' from the specified groups in this
    forcefield. Note that this only removes the specific view
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in this
    view from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::remove(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        if (this->group_remove(mgidx, molview))
        {
            if (not this->group(mgnum).contains(molview.data().number()))
                this->removeFromIndex(mgnum, molview.data().number());
                
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );

        bool removed_mol = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                if (this->group_remove(mgidx, molview))
                {
                    removed_mol = true;
                
                    if (not this->group(mgnum).contains(molview.data().number()))
                        this->removeFromIndex(mgnum, molview.data().number());
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return removed_mol;
    }
}

/** Remove the views of the molecule in 'molviews' from the specified 
    groups in this forcefield. Note that this only removes the specific view
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in these
    views from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::remove(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        ViewsOfMol removed_views = this->group_remove(mgidx, molviews);
        
        if (not removed_views.isEmpty())
        {
            if (not this->group(mgnum).contains(molviews.number()))
                this->removeFromIndex(mgnum, molviews.number());
                
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        bool removed_mol = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                ViewsOfMol removed_views = this->group_remove(mgidx, molviews);
                
                if (not removed_views.isEmpty())
                {
                    removed_mol = true;
                
                    if (not this->group(mgnum).contains(molviews.number()))
                        this->removeFromIndex(mgnum, molviews.number());
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return removed_mol;
    }
}

/** Remove the views of the molecules in 'molecules' from the specified 
    groups in this forcefield. Note that this only removes the specific view
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in these
    views from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::remove(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        QList<ViewsOfMol> removed_mols = this->group_remove(mgidx, molecules);
        QSet<MolNum> removed_molnums;

        const MoleculeGroup &molgroup = this->group(mgnum);

        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            if (not molgroup.contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }
        
        if (not removed_molnums.isEmpty())
        {
            this->removeFromIndex(mgnum, removed_molnums);
            return true;
        }
        else 
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );

        bool mols_removed = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                QList<ViewsOfMol> removed_mols = this->group_remove(mgidx, molecules);
                QSet<MolNum> removed_molnums;
                
                const MoleculeGroup &molgroup = this->group(mgnum);
                
                foreach (const ViewsOfMol &removed_mol, removed_mols)
                {
                    if (not molgroup.contains(removed_mol.number()))
                        removed_molnums.insert(removed_mol.number());
                }
                
                if (not removed_molnums.isEmpty())
                {
                    this->removeFromIndex(mgnum, removed_molnums);
                    mols_removed = true;
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Remove the views of the molecules in the molecule group 'molgroup' f
    from the specified groups in this forcefield. Note that this only 
    removes the specific view (and indeed only the first copy of this view 
    if there are duplicates) - it does not remove the atoms in these
    views from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::remove(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->remove(molgroup.molecules(), mgid);
}

/** Remove the view 'molview' from the specified groups in this
    forcefield. Note that this only removes the specific view
    (and indeed all copies of this view if there 
    are duplicates) - it does not remove the atoms in this
    view from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::removeAll(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        if (this->group_removeAll(mgidx, molview))
        {
            if (not this->group(mgnum).contains(molview.data().number()))
                this->removeFromIndex(mgnum, molview.data().number());
                
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        bool mols_removed = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                if (this->group_removeAll(mgidx, molview))
                {
                    mols_removed = true;
                
                    if (not this->group(mgnum).contains(molview.data().number()))
                        this->removeFromIndex(mgnum, molview.data().number());
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Remove the views of the molecule in 'molviews' from the specified 
    groups in this forcefield. Note that this only removes the specific view
    (and indeed all copies of this view if there 
    are duplicates) - it does not remove the atoms in these
    views from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::removeAll(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        ViewsOfMol removed_views = this->group_removeAll(mgidx, molviews);
        
        if (not removed_views.isEmpty())
        {
            if (not this->group(mgnum).contains(molviews.number()))
                this->removeFromIndex(mgnum, molviews.number());
                
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        bool mols_removed = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                ViewsOfMol removed_views = this->group_removeAll(mgidx, molviews);
                
                if (not removed_views.isEmpty())
                {
                    mols_removed = true;
                
                    if (not this->group(mgnum).contains(molviews.number()))
                        this->removeFromIndex(mgnum, molviews.number());
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Remove the views of the molecules in 'molecules' from the specified 
    groups in this forcefield. Note that this only removes the specific view
    (and indeed all copies of this view if there 
    are duplicates) - it does not remove the atoms in these
    views from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::removeAll(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
    
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        QList<ViewsOfMol> removed_mols = this->group_removeAll(mgidx, molecules);
        QSet<MolNum> removed_molnums;

        const MoleculeGroup &molgroup = this->group(mgnum);

        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            if (not molgroup.contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }
        
        if (not removed_molnums.isEmpty())
        {
            this->removeFromIndex(mgnum, removed_molnums);
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );

        bool mols_removed = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                QList<ViewsOfMol> removed_mols = this->group_removeAll(mgidx, molecules);
                QSet<MolNum> removed_molnums;
                
                const MoleculeGroup &molgroup = this->group(mgnum);
                
                foreach (const ViewsOfMol &removed_mol, removed_mols)
                {
                    if (not molgroup.contains(removed_mol.number()))
                        removed_molnums.insert(removed_mol.number());
                }
                
                if (not removed_molnums.isEmpty())
                {
                    this->removeFromIndex(mgnum, removed_molnums);
                    mols_removed = true;
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Remove the views of the molecules in the molecule group 'molgroup' f
    from the specified groups in this forcefield. Note that this only 
    removes the specific view (and indeed all copies of this view 
    if there are duplicates) - it does not remove the atoms in these
    views from all of the other views
    
    \throw SireMol::missing_group
*/
bool FF::removeAll(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->removeAll(molgroup.molecules(), mgid);
}

/** Completely remove the molecule with number 'molnum' from the 
    forcefield groups identified by 'mgid'
    
    \throw SireMol::missing_group
*/
bool FF::remove(MolNum molnum, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        ViewsOfMol removed_views = this->group_remove(mgidx, molnum);

        if (not removed_views.isEmpty())
        {
            this->removeFromIndex(mgnum, molnum);
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        bool mols_removed = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                ViewsOfMol removed_views = this->group_remove(mgidx, molnum);
                
                if (not removed_views.isEmpty())
                {
                    this->removeFromIndex(mgnum, molnum);
                    mols_removed = true;
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Remove all of the molecules whose numbers are in 'molnums'
    from all of the forcefield groups identified by 'mgid'
    
    \throw SireMol::missing_group
*/
bool FF::remove(const QSet<MolNum> &molnums, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return false;
        
    else if (mgnums.count() == 1)
    {
        MGNum mgnum = *(mgnums.constBegin());
        MGIdx mgidx = this->mgIdx(mgnum);
        
        QList<ViewsOfMol> removed_mols = this->group_remove(mgidx, molnums);
        
        QSet<MolNum> removed_molnums;
        
        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            removed_molnums.insert( removed_mol.number() );
        }
        
        if (not removed_molnums.isEmpty())
        {
            this->removeFromIndex(mgnum, removed_molnums);
            return true;
        }
        else
            return false;
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
        
        bool mols_removed = false;
        
        try
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                
                QList<ViewsOfMol> removed_mols = this->group_remove(mgidx, molnums);
                
                QSet<MolNum> removed_molnums;
                
                foreach (const ViewsOfMol &removed_mol, removed_mols)
                {
                    removed_molnums.insert(removed_mol.number());
                }
                
                if (not removed_molnums.isEmpty())
                {
                    this->removeFromIndex(mgnum, removed_molnums);
                    mols_removed = true;
                }
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Update this forcefield so that it uses the version of the molecule
    that is present in 'moldata' - this does nothing if this molecule
    is not in this forcefield. This uses the existing property names
    to find the updated properties that contain the forcefield parameters
    for this molecule. If you want to change the property names, then you
    must remove the molecule, then re-add it with the new names.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::update(const MoleculeData &moldata, bool auto_commit)
{
    if (not this->contains(moldata.number()))
        return;

    const QList<MGNum> &mgnums = this->groupsContaining(moldata.number());

    BOOST_ASSERT(not mgnums.isEmpty());

    if (auto_commit and this->needsAccepting())
        this->accept();
    
    if (mgnums.count() == 1)
    {
        this->group_update( this->mgIdx(*(mgnums.constBegin())),
                            moldata, auto_commit );
    }
    else
    {
        foreach (MGNum mgnum, mgnums)
        {
            this->group_update( this->mgIdx(mgnum), moldata, auto_commit );
        }
    }
    
    if (auto_commit and this->needsAccepting())
        this->accept();
}

/** Update the data of the molecule that is view in 'molview'. This
    updates all atoms, even those that are not part of the view */
void FF::update(const MoleculeView &molview, bool auto_commit)
{
    this->update(molview.data(), auto_commit);
}

/** Update the molecules in this forcefield so that they have the 
    same version as in 'molecules'. The molecules will use the existing
    property names to find the properties that contain the forcefield
    parameters. If you want to change the property names, then you
    must remove the molecule, then re-add it with the new names
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::update(const Molecules &molecules, bool auto_commit)
{
    if (molecules.isEmpty())
        return;
    
    else if (molecules.count() == 1)
    {
        this->update( molecules.constBegin()->data(), auto_commit );
        return;
    }

    if (auto_commit and this->needsAccepting())
        this->accept();

    //get the numbers of the groups that contain these molecules...
    int ngroups = this->nGroups();
    
    if (ngroups == 1)
    {
        this->group_update(0, molecules, auto_commit);
    }
    else
    {
        QSet<MGNum> mgnums;

        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (this->contains(it->number()))
            {
                const QList<MGNum> &mol_mgnums = this->groupsContaining(it->number());
            
                foreach (MGNum mgnum, mol_mgnums)
                {
                    mgnums.insert(mgnum);
                }
                
                if (mgnums.count() == ngroups)
                    break;
            }
        }
        
        if (mgnums.count() == 1)
        {
            //only one group needs to be updated
            MGNum mgnum = *(mgnums.constBegin());
            MGIdx mgidx = this->mgIdx(mgnum);
            
            this->group_update(mgidx, molecules, auto_commit);
        }
        else
        {
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                this->group_update(mgidx, molecules, auto_commit);
            }
        }
    }
    
    if (auto_commit and this->needsAccepting())
        this->accept();
}

/** Update this forcefield so that it has the same version molecules
    as those contained in 'molgroup'. 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::update(const MoleculeGroup &molgroup, bool auto_commit)
{
    this->update(molgroup.molecules(), auto_commit);
}

/** Set the contents of the forcefield groups identified by 'mgid'
    so that they only contain the molecule viewed in 'molview'.
    This will use the version of the molecule that exists already
    in this forcefield, not the version in 'molview'
    
    This will use the supplied map to find the property names
    of the parameters required by this forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const MoleculeView &molview, 
                     const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;

    if (this->nGroups() == 1)
    {
        //there is only one group and it is being updated
        MGNum mgnum = *(mgnums.constBegin());
        
        if (this->group_setContents(0, molview, map))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, molview.data().number());
        }
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
    
        try
        {
            //first clear any existing molecules
            foreach (MGNum mgnum, mgnums)
            {
                this->group_removeAll(this->mgIdx(mgnum));
                this->clearIndex(mgnum);
            }
            
            //now get the current version of the molecule
            PartialMolecule view(molview);
            view.update( this->matchToExistingVersion(view.data()) );
            
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                this->group_setContents(mgidx, view, map);
                
                this->addToIndex(mgnum, view.number());
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}
                     
/** Set the contents of the forcefield groups identified by 'mgid'
    so that they only contain the molecule views in 'molviews'.
    This will use the version of the molecule that exists already
    in this forcefield, not the version in 'molview'
    
    This will use the supplied map to find the property names
    of the parameters required by this forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const ViewsOfMol &molviews, 
                     const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;

    if (this->nGroups() == 1)
    {
        //there is only one group and it is being updated
        MGNum mgnum = *(mgnums.constBegin());
        
        if (this->group_setContents(0, molviews, map))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, molviews.number());
        }
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
    
        try
        {
            //first clear any existing molecules
            foreach (MGNum mgnum, mgnums)
            {
                this->group_removeAll(this->mgIdx(mgnum));
                this->clearIndex(mgnum);
            }
            
            //now get the current version of the molecule
            ViewsOfMol views(molviews);
            views.update( this->matchToExistingVersion(views.data()) );
            
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                this->group_setContents(mgidx, views, map);
                
                this->addToIndex(mgnum, views.number());
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}

/** Set the contents of the forcefield groups identified by 'mgid'
    so that they only contain the molecules in 'molecules'.
    This will use the version of the molecule that exists already
    in this forcefield, not the version in 'molview'
    
    This will use the supplied map to find the property names
    of the parameters required by this forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const Molecules &molecules, 
                     const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.isEmpty())
        return;

    if (this->nGroups() == 1)
    {
        //there is only one group and it is being updated
        MGNum mgnum = *(mgnums.constBegin());
        
        if (this->group_setContents(0, molecules, map))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, molecules.molNums());
        }
    }
    else
    {
        boost::shared_ptr<FF> old_state( this->clone() );
    
        try
        {
            //first clear any existing molecules
            foreach (MGNum mgnum, mgnums)
            {
                this->group_removeAll(this->mgIdx(mgnum));
                this->clearIndex(mgnum);
            }
            
            //now get the current version of the molecules
            Molecules mols = this->matchToExistingVersion(molecules);
            
            QSet<MolNum> molnums = mols.molNums();
            
            foreach (MGNum mgnum, mgnums)
            {
                MGIdx mgidx = this->mgIdx(mgnum);
                this->group_setContents(mgidx, mols, map);
                
                this->addToIndex(mgnum, molnums);
            }
        }
        catch(...)
        {
            this->copy(*old_state);
            throw;
        }
    }
}

/** Set the contents of the forcefield groups identified by 'mgid'
    so that they only contain the molecules contained in 'molgroup'.
    This will use the version of the molecule that exists already
    in this forcefield, not the version in 'molview'
    
    This will use the supplied map to find the property names
    of the parameters required by this forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const MoleculeGroup &molgroup, 
                     const PropertyMap &map)
{
    this->setContents(mgid, molgroup.molecules(), map);
}

/** Add the passed view of the molecule to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const MoleculeView &molview, const MGID &mgid)
{
    this->add(molview, mgid, PropertyMap());
}

/** Add the passed views of the molecule to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->add(molviews, mgid, PropertyMap());
}

/** Add the passed molecules to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const Molecules &molecules, const MGID &mgid)
{
    this->add(molecules, mgid, PropertyMap());
}

/** Add the molecules in the passed MoleculeGroup to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::add(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->add(molgroup, mgid, PropertyMap());
}

/** Add the passed view of the molecule to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    Only add this view to groups that don't already contain
    this view (the whole view, not part of it)
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const MoleculeView &molview, const MGID &mgid)
{
    this->addIfUnique(molview, mgid, PropertyMap());
}

/** Add the passed views of the molecule to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    Only add views to groups that don't already contain
    them (the whole view, not part of it, and can add some views)
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->addIfUnique(molviews, mgid, PropertyMap());
}

/** Add the passed molecules to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    Only add the views of molecules to groups that don't already contain
    them (the whole view, not part of it)
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const Molecules &molecules, const MGID &mgid)
{
    this->addIfUnique(molecules, mgid, PropertyMap());
}

/** Add the molecules in the passed MoleculeGroup to the molecule groups 
    identified by 'mgid' using the default properties to
    find the parameters needed by this forcefield
    
    Only add the views of molecules to groups that don't already contain
    them (the whole view, not part of it)
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->addIfUnique(molgroup, mgid, PropertyMap());
}

/** Set the contents of this forcefield to just contain 'molview', 
    using the default locations to find the properties that contain
    the forcefield parameters for this molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const MoleculeView &molview)
{
    this->setContents(mgid, molview, PropertyMap());
}

/** Set the contents of this forcefield to the molecule views in 'molviews', 
    using the default locations to find the properties that contain
    the forcefield parameters for this molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const ViewsOfMol &molviews)
{
    this->setContents(mgid, molviews, PropertyMap());
}

/** Set the contents of this forcefield to the molecules 'molecules', 
    using the default locations to find the properties that contain
    the forcefield parameters for this molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const Molecules &molecules)
{
    this->setContents(mgid, molecules, PropertyMap());
}
                         
/** Set the contents of this forcefield to contains the molecules in 'molgroup', 
    using the default locations to find the properties that contain
    the forcefield parameters for this molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FF::setContents(const MGID &mgid, const MoleculeGroup &molgroup)
{
    this->setContents(mgid, molgroup, PropertyMap());
}

/** Return the names of all of the properties available to this forcefield */
QStringList FF::propertyKeys() const
{
    return this->properties().propertyKeys();
}
