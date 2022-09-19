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

#include <QMutex>

#include "moleculegroups.h"

#include "molecules.h"
#include "molecule.h"
#include "viewsofmol.h"
#include "partialmolecule.h"
#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "mover.hpp"
#include "editor.hpp"

#include "mgnum.h"
#include "mgidx.h"
#include "mgname.h"

#include "molnum.h"
#include "molidx.h"
#include "molname.h"

#include "segid.h"
#include "chainid.h"
#include "resid.h"
#include "cgid.h"
#include "atomid.h"

#include "select.h"
#include "selector.hpp"

#include "tostring.h"

#include "SireBase/slice.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

/////////////
///////////// Implementation of MolGroupsBase
/////////////

static const RegisterMetaType<MolGroupsBase> r_molgroupsbase(MAGIC_ONLY,
                                                "SireMol::MolGroupsBase");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const MolGroupsBase &molgroupsbase)
{
    writeHeader(ds, r_molgroupsbase, 2);

    SharedDataStream sds(ds);

    sds << molgroupsbase.mgidx_to_num << static_cast<const Property&>(molgroupsbase);

    return ds;
}

/** Read from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       MolGroupsBase &molgroupsbase)
{
    VersionID v = readHeader(ds, r_molgroupsbase);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> molgroupsbase.mgidx_to_num >> static_cast<Property&>(molgroupsbase);

        molgroupsbase.reindex();
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> molgroupsbase.mgidx_to_num
            >> molgroupsbase.mgname_to_mgnum
            >> molgroupsbase.molnum_to_mgnum
            >> static_cast<Property&>(molgroupsbase);
    }
    else
        throw version_error(v, "1", r_molgroupsbase, CODELOC);

    return ds;
}

/** Null constructor */
MolGroupsBase::MolGroupsBase() : Property()
{}

/** Copy constructor */
MolGroupsBase::MolGroupsBase(const MolGroupsBase &other)
              : Property(other),
                mgidx_to_num(other.mgidx_to_num),
                mgname_to_mgnum(other.mgname_to_mgnum),
                molnum_to_mgnum(other.molnum_to_mgnum)
{}

/** Destructor */
MolGroupsBase::~MolGroupsBase()
{}

/** Copy assignment operator */
MolGroupsBase& MolGroupsBase::operator=(const MolGroupsBase &other)
{
    if (this != &other)
    {
        mgidx_to_num = other.mgidx_to_num;
        mgname_to_mgnum = other.mgname_to_mgnum;
        molnum_to_mgnum = other.molnum_to_mgnum;
        Property::operator=(other);
    }

    return *this;
}

/** Return a const reference to the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& MolGroupsBase::operator[](MGNum mgnum) const
{
    return this->at(mgnum);
}

/** Return a const reference to the molecule group with name 'mgname'

    \throw SireMol::missing_group
*/
const MoleculeGroup& MolGroupsBase::operator[](const MGName &mgname) const
{
    return this->at(mgname);
}

/** Return a const reference to the molecule group at index 'mgidx'

    \throw SireError::invalid_index
*/
const MoleculeGroup& MolGroupsBase::operator[](MGIdx mgidx) const
{
    return this->at(mgidx);
}

/** Return a const reference to the molecule group that is
    identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
const MoleculeGroup& MolGroupsBase::operator[](const MGID &mgid) const
{
    return this->at(mgid);
}

/** Return all of the views of the molecule with number 'molnum'
    that are contained in this set of groups. Note that if the
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views

    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::operator[](MolNum molnum) const
{
    return this->at(molnum);
}

ViewsOfMol MolGroupsBase::operator[](int i) const
{
    return this->at(MolIdx(i));
}

ViewsOfMol MolGroupsBase::operator[](const QString &name) const
{
    return this->at(MolName(name));
}

QList<MolViewPtr> MolGroupsBase::operator[](const SireBase::Slice &slice) const
{
    const auto molnums = this->molNums();

    QList<MolViewPtr> views;

    for (auto it = slice.begin(molnums.count()); not it.atEnd(); it.next())
    {
        views.append(this->operator[](molnums.at(it.value())));
    }

    return views;
}

/** Return all of the views of the molecule identified by 'molid'
    that are contained in this set of groups. Note that if the
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol MolGroupsBase::operator[](const MolID &molid) const
{
    return this->at(molid);
}

Segment MolGroupsBase::operator[](const SegID &segid) const
{
    return this->at(segid);
}

Chain MolGroupsBase::operator[](const ChainID &chainid) const
{
    return this->at(chainid);
}

Residue MolGroupsBase::operator[](const ResID &resid) const
{
    return this->at(resid);
}

CutGroup MolGroupsBase::operator[](const CGID &cgid) const
{
    return this->at(cgid);
}

Atom MolGroupsBase::operator[](const AtomID &atomid) const
{
    return this->at(atomid);
}

/** Get the number of the molecule group whose number is 'mgnum'.
    This is an obvious function, only provided as a shortcut
    to prevent the MGID function being called if an MGNum is passed.

    \throw SireMol::missing_group
*/
MGNum MolGroupsBase::getGroupNumber(MGNum mgnum) const
{
    this->assertContains(mgnum);
    return mgnum;
}

/** Return the number of the molecule group that is called 'mgname'.

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
MGNum MolGroupsBase::getGroupNumber(const MGName &mgname) const
{
    QList<MGNum> mgnums = this->map(mgname);

    if (mgnums.count() > 1)
        throw SireMol::duplicate_group( QObject::tr(
            "There are multiple molecule groups called \"%1\" "
            "in this set - with group numbers %2.")
                .arg(mgname).arg(Sire::toString(mgnums)), CODELOC );

    return mgnums.first();
}

/** Return the number of the group at index 'mgidx'

    \throw SireError::invalid_index
*/
MGNum MolGroupsBase::getGroupNumber(MGIdx mgidx) const
{
    return mgidx_to_num.at( mgidx.map(mgidx_to_num.count()) );
}

/** Return the number of the groups that matches 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
MGNum MolGroupsBase::getGroupNumber(const MGID &mgid) const
{
    QList<MGNum> mgnums = this->map(mgid);

    if (mgnums.count() > 1)
        throw SireMol::duplicate_group( QObject::tr(
            "There is more than one molecule group that matches "
            "the ID \"%1\". Matching groups have numbers %2.")
                .arg(mgid.toString()).arg(Sire::toString(mgnums)),
                    CODELOC );

    return mgnums.first();
}

/** Return the index of the group with number 'mgnum'

    \throw SireMol::missing_group
*/
MGIdx MolGroupsBase::mgIdx(MGNum mgnum) const
{
    int i = mgidx_to_num.indexOf(mgnum);

    if (i == -1)
        throw SireMol::missing_group( QObject::tr(
            "There is no molecule group with number %1 in this set. "
            "Available groups are %2.")
                .arg(mgnum).arg( Sire::toString(mgidx_to_num) ), CODELOC );

    return MGIdx(i);
}

/** Return the numbers of all groups in this set that are called
    'mgname'

    \throw SireMol::missing_group
*/
QList<MGNum> MolGroupsBase::map(const MGName &mgname) const
{
    QList<MGNum> mgnums;

    if (mgname.isCaseSensitive())
    {
        QHash< MGName,QList<MGNum> >::const_iterator it = mgname_to_mgnum.find(mgname);

        if (it != mgname_to_mgnum.end())
            mgnums = it.value();
    }
    else
    {
        QString lower_name = QString(mgname).toLower();

        for (QHash< MGName,QList<MGNum> >::const_iterator
                                                it = mgname_to_mgnum.constBegin();
             it != mgname_to_mgnum.constEnd();
             ++it)
        {
            if (QString(it.key()).toLower() == lower_name)
                mgnums += it.value();
        }
    }

    if (mgnums.isEmpty())
        throw SireMol::missing_group( QObject::tr(
            "There are no molecule groups called \"%1\" in this set. "
            "Available groups are %2.")
                .arg(mgname).arg(Sire::toString(mgname_to_mgnum.keys())),
                    CODELOC );

    return mgnums;
}

/** Return the list of numbers of groups that have the number 'mgnum'.
    This is a simple and obvious function that acts as a shortcut
    preventing map(const MGID&) being called for an MGNum

    \throw SireMol::missing_group
*/
QList<MGNum> MolGroupsBase::map(MGNum mgnum) const
{
    this->assertContains(mgnum);

    QList<MGNum> mgnums;
    mgnums.append(mgnum);

    return mgnums;
}

/** Return the list (of only one) molecule group that is at
    index 'mgidx'

    \throw SireError::invalid_index
*/
QList<MGNum> MolGroupsBase::map(MGIdx mgidx) const
{
    QList<MGNum> mgnums;

    mgnums.append( this->getGroupNumber(mgidx) );

    return mgnums;
}

/** Map the molecule group ID 'mgid' to the list of molecule
    group numbers of the groups that match this ID in this set.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
QList<MGNum> MolGroupsBase::map(const MGID &mgid) const
{
    return mgid.map(*this);
}

/** Simple function that just checks if a molecule with number
    'molnum' is in the set, and returns it. This shortcuts
    the getMoleculeNumber(const MolID&) function in the case
    of MolNums

    \throw SireMol::missing_molecule
*/
MolNum MolGroupsBase::getMoleculeNumber(MolNum molnum) const
{
    this->assertContains(molnum);
    return molnum;
}

/** Return the number of the molecule at index 'molidx' in
    this set

    \throw SireError::invalid_index
*/
MolNum MolGroupsBase::getMoleculeNumber(MolIdx molidx) const
{
    //map this to a valid index
    int n = molidx.map( molnum_to_mgnum.count() );

    const auto molgroups = this->groups();

    if (molgroups.count() == 1 or n < molgroups[0].read().nMolecules())
    {
        return molgroups[0].read().at(molidx).data().number();
    }
    else
    {
        //we now run through each molecule of each group, discounting
        //any that have been seen previously
        n -= molgroups[0].read().nMolecules();

        QSet<MGNum> seen_groups;
        seen_groups.insert(molgroups[0].read().number());

        for (int i=1; i<molgroups.count(); ++i)
        {
            const auto molgroup = molgroups[i].read();

            for (int j=0; j<molgroup.nMolecules(); ++j)
            {
                MolNum molnum = molgroup[MolIdx(j)].data().number();

                bool seen_before = false;
                for (auto mgnum : molnum_to_mgnum.value(molnum))
                {
                    if (seen_groups.contains(mgnum))
                    {
                        seen_before = true;
                        break;
                    }
                }

                if (not seen_before)
                {
                    n -= 1;

                    if (n == 0)
                    {
                        //we have found the molecule
                        return molnum;
                    }
                }
            }

            seen_groups.insert(molgroup.number());
        }
    }

    throw SireError::program_bug( QObject::tr(
            "Could not find the molecule at index %1 despite there being "
            "%2 molecules in the set of molecule groups???")
                .arg(molidx.value()).arg(molnum_to_mgnum.count()), CODELOC );

    return MolNum();
}

/** Return the list of molecule numbers in molidx order */
QList<MolNum> MolGroupsBase::getMoleculeNumbers() const
{
    QList<MolNum> molnums;

    const auto molgroups = this->groups();

    if (molgroups.isEmpty())
        return molnums;

    //first add the molecule numbers from the first group
    molnums = molgroups[0].read().molNums().toList();

    //we now run through each molecule of each group, discounting
    //any that have been seen previously
    QSet<MGNum> seen_groups;
    seen_groups.insert(molgroups[0].read().number());

    for (int i=1; i<molgroups.count(); ++i)
    {
        const auto molgroup = molgroups[i].read();

        for (int j=0; j<molgroup.nMolecules(); ++j)
        {
            MolNum molnum = molgroup[MolIdx(j)].data().number();

            bool seen_before = false;
            for (auto mgnum : molnum_to_mgnum.value(molnum))
            {
                if (seen_groups.contains(mgnum))
                {
                    seen_before = true;
                    break;
                }
            }

            if (not seen_before)
                molnums.append(molnum);
        }

        seen_groups.insert(molgroup.number());
    }

    return molnums;
}

/** Return the list of molecule numbers in molidx order */
QList<MolNum> MolGroupsBase::molNums() const
{
    return this->getMoleculeNumbers();
}

/** Return the number of the molecule called 'molname' from this set.

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
MolNum MolGroupsBase::getMoleculeNumber(const MolName &molname) const
{
    QList<MolNum> molnums = this->map(molname);

    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule with the name \"%1\" "
            "in this set of molecule groups. Molecules with this "
            "name have numbers %2.")
                .arg(molname).arg(Sire::toString(molnums)), CODELOC );

    return molnums.first();
}

/** Return the number of the molecule that matches the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
MolNum MolGroupsBase::getMoleculeNumber(const MolID &molid) const
{
    QList<MolNum> molnums = this->map(molid);

    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule that matches the "
            "ID \"%1\". Molecules that match have number %2.")
                .arg(molid.toString()).arg(Sire::toString(molnums)),
                    CODELOC );

    return molnums.first();
}

/** Return the version number of the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
quint64 MolGroupsBase::getMoleculeVersion(MolNum molnum) const
{
    QHash< MolNum, QList<MGNum> >::const_iterator it = molnum_to_mgnum.constFind(molnum);

    if (it == molnum_to_mgnum.constEnd())
        throw SireMol::missing_molecule( QObject::tr(
                "There is no molecule with number %1 available in this set of groups.")
                    .arg(molnum.toString()), CODELOC );

    return this->getGroup( it.value().at(0) ).getMoleculeVersion(molnum);
}

/** Return the version number of the molecule with ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
quint64 MolGroupsBase::getMoleculeVersion(const MolID &molid) const
{
    return this->getMoleculeVersion( this->getMoleculeNumber(molid) );
}

/** Simple function that provides a shortcut for map(const MolID&)

    \throw SireMol::missing_molecule
*/
QList<MolNum> MolGroupsBase::map(MolNum molnum) const
{
    this->assertContains(molnum);

    QList<MolNum> molnums;
    molnums.append(molnum);

    return molnums;
}

/** Return the number of the molecule at index 'molidx'

    \throw SireError::invalid_index
*/
QList<MolNum> MolGroupsBase::map(MolIdx molidx) const
{
    QList<MolNum> molnums;
    molnums.append( this->getMoleculeNumber(molidx) );

    return molnums;
}

/** Return the numbers of all of the molecules that have the
    name 'molname'

    \throw SireMol::missing_molecule
*/
QList<MolNum> MolGroupsBase::map(const MolName &molname) const
{
    QList<MolNum> molnums;

    //loop over all of the groups in this set
    const QHash<MGNum,const MoleculeGroup*> molgroups = this->getGroups();

    for (QHash<MGNum,const MoleculeGroup*>::const_iterator it = molgroups.constBegin();
         it != molgroups.constEnd();
         ++it)
    {
        try
        {
            molnums += (*it)->map(molname);
        }
        catch(...)
        {}
    }

    if (molnums.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "There are no molecules called \"%1\" in any of the "
            "molecule groups in this set.")
                .arg(molname), CODELOC );

    else if (molnums.count() > 1)
    {
        //remove any duplicates
        molnums = convert_to_qset(molnums).values();
    }

    return molnums;
}

/** Return the numbers of all molecules that match the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> MolGroupsBase::map(const MolID &molid) const
{
    return molid.map(*this);
}

/** Return a const reference to the group at index 'mgidx'

    \throw SireError::invalid_index
*/
const MoleculeGroup& MolGroupsBase::at(MGIdx mgidx) const
{
    return this->at( this->getGroupNumber(mgidx) );
}

/** Return a const reference to the group in this set
    called 'mgname'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
const MoleculeGroup& MolGroupsBase::at(const MGName &mgname) const
{
    return this->at( this->getGroupNumber(mgname) );
}

/** Return a const reference to the group in this set that
    is identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const MoleculeGroup& MolGroupsBase::at(const MGID &mgid) const
{
    return this->at( this->getGroupNumber(mgid) );
}

/** Return all of the views of the molecule with number 'molnum'
    that are contained in this set of groups. Note that if the
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views

    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::at(MolNum molnum) const
{
    //get the list of groups that contain this molecule
    const QList<MGNum> &mgnums = this->groupsContaining(molnum);

    if (mgnums.count() == 1)
        return this->at( mgnums.first() ).molecule(molnum);
    else
    {
        ViewsOfMol mol = this->at( mgnums.first() ).molecule(molnum);

        for (int i=1; i<mgnums.count(); ++i)
        {
            mol.add( this->at( mgnums.at(i) ).molecule(molnum) );
        }

        return mol;
    }
}

/** Return all of the views of the molecule identified by 'molid'
    that are contained in this set of groups. Note that if the
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol MolGroupsBase::at(const MolID &molid) const
{
    QList<MolNum> molnums = molid.map(*this);

    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule that matches the ID "
            "\"%1\". Matching molecules have numbers %2.")
                .arg(molid.toString())
                .arg(Sire::toString(molnums)), CODELOC );

    return this->at(molnums.first());
}

/** Return the segment from this set that matches the ID 'segid'.
    This segment must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
Segment MolGroupsBase::at(const SegID &segid) const
{
    return segid.selectFrom(*this);
}

/** Return the chain from this set that matches the ID 'chainid'.
    This chain must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
Chain MolGroupsBase::at(const ChainID &chainid) const
{
    return chainid.selectFrom(*this);
}

/** Return the residue from this set that matches the ID 'resid'.
    This residue must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
Residue MolGroupsBase::at(const ResID &resid) const
{
    return resid.selectFrom(*this);
}

/** Return the CutGroup from this set that matches the ID 'cgid'.
    This CutGroup must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CutGroup MolGroupsBase::at(const CGID &cgid) const
{
    return cgid.selectFrom(*this);
}

/** Return the atom from this set that matches the ID 'atomid'.
    This atom must be contained in one of the groups in this set.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom MolGroupsBase::at(const AtomID &atomid) const
{
    return atomid.selectFrom(*this);
}

/** Return the MoleculeGroup that matches the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const MoleculeGroup& MolGroupsBase::select(const MGID &mgid) const
{
    return this->at(mgid);
}

/** Return all of the views of the molecule with number 'molnum'
    that are contained in this set of groups. Note that if the
    same view appears in multiple groups, then it will be returned
    multiple times in the returned set of views

    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::select(const MolID &molid) const
{
    return this->at(molid);
}

/** Return the segment from this set that matches the ID 'segid'.
    This segment must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
Segment MolGroupsBase::select(const SegID &segid) const
{
    return this->at(segid);
}

/** Return the chain from this set that matches the ID 'chainid'.
    This chain must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
Chain MolGroupsBase::select(const ChainID &chainid) const
{
    return this->at(chainid);
}

/** Return the residue from this set that matches the ID 'resid'.
    This residue must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
Residue MolGroupsBase::select(const ResID &resid) const
{
    return this->at(resid);
}

/** Return the CutGroup from this set that matches the ID 'cgid'.
    This CutGroup must be wholly contained by one of the groups
    in this set

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CutGroup MolGroupsBase::select(const CGID &cgid) const
{
    return this->at(cgid);
}

/** Return the atom from this set that matches the ID 'atomid'.
    This atom must be contained in one of the groups in this set.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom MolGroupsBase::select(const AtomID &atomid) const
{
    return this->at(atomid);
}

/** Return the result of searching this object for 'search_string' */
SelectResult MolGroupsBase::search(const QString &search_string) const
{
    return Select(search_string)(*this);
}

/** Return a list of all of the molecule groups in this set */
QList<MolGroupPtr> MolGroupsBase::selectAll() const
{
    QList<MolGroupPtr> molgroups;

    QHash<MGNum,const MoleculeGroup*> groups = this->getGroups();

    foreach (MGNum mgnum, mgidx_to_num)
    {
        molgroups.append( *(groups.value(mgnum)) );
    }

    return molgroups;
}

/** Return a list of all of the molecule groups in this set */
QList<MolGroupPtr> MolGroupsBase::groups() const
{
    return this->selectAll();
}

/** Return a list of the numbers of all of the groups in this set */
QList<MGNum> MolGroupsBase::groupNumbers() const
{
    return mgidx_to_num;
}

/** Return a list of the names of all of the groups in this set */
QList<MGName> MolGroupsBase::groupNames() const
{
    return this->mgNames();
}

/** Obvious shortcut for select(const MGID&)

    \throw SireMol::missing_group
*/
QList<MolGroupPtr> MolGroupsBase::selectAll(MGNum mgnum) const
{
    QList<MolGroupPtr> molgroups;
    molgroups.append( this->at(mgnum) );

    return molgroups;
}

/** Obvious shortcut for select(const MGID&)

    \throw SireError::invalid_index
*/
QList<MolGroupPtr> MolGroupsBase::selectAll(MGIdx mgidx) const
{
    QList<MolGroupPtr> molgroups;
    molgroups.append( this->at(mgidx) );

    return molgroups;
}

/** Return all of the molecule groups that are called 'mgname'

    \throw SireMol::missing_group
*/
QList<MolGroupPtr> MolGroupsBase::selectAll(const MGName &mgname) const
{
    QHash< MGName,QList<MGNum> >::const_iterator it = mgname_to_mgnum.find(mgname);

    if (it == mgname_to_mgnum.end())
        throw SireMol::missing_group( QObject::tr(
            "There are no groups in this set called \"%1\". "
            "The groups in this set are called %2.")
                .arg(mgname).arg(Sire::toString(mgname_to_mgnum.keys())),
                    CODELOC );

    //now get all of the groups
    QVarLengthArray<const MoleculeGroup*,10> groups;
    this->getGroups(*it, groups);

    QList<MolGroupPtr> molgroups;

    int ngroups = groups.count();
    const MoleculeGroup* const *groups_array = groups.constData();

    for (int i=0; i<ngroups; ++i)
    {
        molgroups.append( MoleculeGroup( *(groups_array[i]) ) );
    }

    return molgroups;
}

/** Return all of the molecule groups that match the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
QList<MolGroupPtr> MolGroupsBase::selectAll(const MGID &mgid) const
{
    //get the list of numbers that match this ID
    QList<MGNum> mgnums = this->map(mgid);

    //now pick up those groups...
    QVarLengthArray<const MoleculeGroup*,10> groups;
    this->getGroups(mgnums, groups);

    QList<MolGroupPtr> molgroups;

    int ngroups = groups.count();
    const MoleculeGroup* const *groups_array = groups.constData();

    for (int i=0; i<ngroups; ++i)
    {
        molgroups.append( MoleculeGroup( *(groups_array[i]) ) );
    }

    return molgroups;
}

/** Return the views of the molecule(s) that match the molecule ID
    'molid'. This returns all views of the molecule in the groups,
    and if a view is contained multiple times, then multiple copies
    of that view will be returned.

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<ViewsOfMol> MolGroupsBase::selectAll(const MolID &molid) const
{
    //get the numbers of molecules that match this ID
    QList<MolNum> molnums = this->map(molid);

    QList<ViewsOfMol> molviews;

    foreach (MolNum molnum, molnums)
    {
        molviews.append(this->at(molnum));
    }

    return molviews;
}

/** Return all of the segments from this set that match the ID 'segid'.
    The returned segments are arranged by molecule, and only one copy
    of each segment is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Segment> > MolGroupsBase::selectAll(const SegID &segid) const
{
    return segid.selectAllFrom(*this);
}

/** Return all of the chains from this set that match the ID 'chainid'.
    The returned chains are arranged by molecule, and only one copy
    of each chain is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Chain> > MolGroupsBase::selectAll(const ChainID &chainid) const
{
    return chainid.selectAllFrom(*this);
}

/** Return all of the residues from this set that match the ID 'resid'.
    The returned residues are arranged by molecule, and only one copy
    of each residue is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Residue> > MolGroupsBase::selectAll(const ResID &resid) const
{
    return resid.selectAllFrom(*this);
}

/** Return all of the CutGroups from this set that match the ID 'cgid'.
    The returned CutGroups are arranged by molecule, and only one copy
    of each CutGroup is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<CutGroup> > MolGroupsBase::selectAll(const CGID &cgid) const
{
    return cgid.selectAllFrom(*this);
}

/** Return all of the atoms from this set that match the ID 'atomid'.
    The returned atoms are arranged by molecule, and only one copy
    of each atom is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Atom> >
MolGroupsBase::selectAll(const AtomID &atomid) const
{
    return atomid.selectAllFrom(*this);
}

/** Return the molecule group that has number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& MolGroupsBase::group(MGNum mgnum) const
{
    return this->at(mgnum);
}

/** Return the molecule group that has name 'mgname'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
*/
const MoleculeGroup& MolGroupsBase::group(const MGName &mgname) const
{
    return this->at(mgname);
}

/** Return the molecule group at index 'mgidx'

    \throw SireError::invalid_index
*/
const MoleculeGroup& MolGroupsBase::group(MGIdx mgidx) const
{
    return this->at(mgidx);
}

/** Return the molecule group that matches the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const MoleculeGroup& MolGroupsBase::group(const MGID &mgid) const
{
    return this->at(mgid);
}

/** Obvious shortcut for groups(const MGID&)

    \throw SireMol::missing_group
*/
QList<MolGroupPtr> MolGroupsBase::groups(MGNum mgnum) const
{
    return this->selectAll(mgnum);
}

/** Obvious shortcut for groups(const MGID&)

    \throw SireMol::invalid_index
*/
QList<MolGroupPtr> MolGroupsBase::groups(MGIdx mgidx) const
{
    return this->selectAll(mgidx);
}

/** Return all of the groups called 'mgname'

    \throw SireMol::missing_group
*/
QList<MolGroupPtr> MolGroupsBase::groups(const MGName &mgname) const
{
    return this->selectAll(mgname);
}

/** Return all of the groups that match the ID 'mgid'

    \throw SireMol::missing_group
*/
QList<MolGroupPtr> MolGroupsBase::groups(const MGID &mgid) const
{
    return this->selectAll(mgid);
}

/** Return all of the views of the molecule that has number 'molnum'

    \throw SireMol::missing_molecule
*/
ViewsOfMol MolGroupsBase::molecule(MolNum molnum) const
{
    return this->at(molnum);
}

/** Return all of the views of the molecule that matches 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol MolGroupsBase::molecule(const MolID &molid) const
{
    return this->at(molid);
}

/** Obvious shortcut for molecules(const MolID&)

    \throw SireMol::missing_molecule
*/
QList<ViewsOfMol> MolGroupsBase::molecules(MolNum molnum) const
{
    return this->selectAll(molnum);
}

/** Return all of the molecules that match the ID 'molid'

    \throw SireMol::missing_molecule
*/
QList<ViewsOfMol> MolGroupsBase::molecules(const MolID &molid) const
{
    return this->selectAll(molid);
}

/** Return the segment that matches the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
Segment MolGroupsBase::segment(const SegID &segid) const
{
    return this->at(segid);
}

/** Return the chain that matches the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
Chain MolGroupsBase::chain(const ChainID &chainid) const
{
    return this->at(chainid);
}

/** Return the residue that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
Residue MolGroupsBase::residue(const ResID &resid) const
{
    return this->at(resid);
}

/** Return the CutGroup that matches the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CutGroup MolGroupsBase::cutGroup(const CGID &cgid) const
{
    return this->at(cgid);
}

/** Return the atom that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom MolGroupsBase::atom(const AtomID &atomid) const
{
    return this->at(atomid);
}

/** Return all of the segments from this set that match the ID 'segid'.
    The returned segments are arranged by molecule, and only one copy
    of each segment is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Segment> > MolGroupsBase::segments(const SegID &segid) const
{
    return this->selectAll(segid);
}

/** Return all of the chains from this set that match the ID 'chainid'.
    The returned chains are arranged by molecule, and only one copy
    of each chain is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Chain> > MolGroupsBase::chains(const ChainID &chainid) const
{
    return this->selectAll(chainid);
}

/** Return all of the residues from this set that match the ID 'resid'.
    The returned residues are arranged by molecule, and only one copy
    of each residue is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Residue> > MolGroupsBase::residues(const ResID &resid) const
{
    return this->selectAll(resid);
}

/** Return all of the CutGroups from this set that match the ID 'cgid'.
    The returned CutGroups are arranged by molecule, and only one copy
    of each CutGroup is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<CutGroup> > MolGroupsBase::cutGroups(const CGID &cgid) const
{
    return this->selectAll(cgid);
}

/** Return all of the atoms from this set that match the ID 'atomid'.
    The returned atoms are arranged by molecule, and only one copy
    of each atom is returned, regardless of how many times it appears
    in this set.

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QHash< MolNum,Selector<Atom> > MolGroupsBase::atoms(const AtomID &atomid) const
{
    return this->selectAll(atomid);
}

/** Return whether or not this set contains the group with number 'mgnum' */
bool MolGroupsBase::contains(MGNum mgnum) const
{
    return mgidx_to_num.contains(mgnum);
}

/** Return whether any of the groups contain any view of the molecule
    with number 'molnum' */
bool MolGroupsBase::contains(MolNum molnum) const
{
    return molnum_to_mgnum.contains(molnum);
}

/** Return whether any of the groups contains any of the molecules whose
    numbers are in 'molnums' */
bool MolGroupsBase::contains(const QList<MolNum> &molnums) const
{
    if (molnum_to_mgnum.isEmpty())
        return false;

    for (QList<MolNum>::const_iterator it = molnums.constBegin();
         it != molnums.constEnd();
         ++it)
    {
        if (molnum_to_mgnum.contains(*it))
            return true;
    }

    return false;
}

/** Return whether or not any of the groups contains the view 'molview' */
bool MolGroupsBase::contains(const MoleculeView &molview) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator
                                  it = molnum_to_mgnum.find(molview.data().number());

    if (it == molnum_to_mgnum.end())
        return false;

    QVarLengthArray<const MoleculeGroup*,10> groups;

    this->getGroups(*it, groups);

    int count = groups.count();
    const MoleculeGroup* const *groups_array = groups.constData();

    for (int i=0; i<count; ++i)
    {
        if (groups_array[i]->contains(molview))
            return true;
    }

    return false;
}

/** Return whether or not this set contains all of the views of
    the molecule in 'molviews'. The views can be contained in
    multiple groups. */
bool MolGroupsBase::contains(const ViewsOfMol &molviews) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator
                                    it = molnum_to_mgnum.find(molviews.number());

    if (it == molnum_to_mgnum.end())
        return false;

    QVarLengthArray<const MoleculeGroup*,10> groups;

    this->getGroups(*it, groups);

    for (int i=0; i<molviews.nViews(); ++i)
    {
        PartialMolecule view = molviews.valueAt(i);

        bool found_view = false;

        int count = groups.count();
        const MoleculeGroup* const *groups_array = groups.constData();

        for (int i=0; i<count; ++i)
        {
            if (groups_array[i]->contains(view))
            {
                found_view = true;
                break;
            }
        }

        if (not found_view)
            return false;
    }

    return true;
}

/** Return whether or not this set of groups contains all of the views
    of all of the molecules in 'molecules'. These views can be spread
    over lots of groups */
bool MolGroupsBase::contains(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        if (not this->contains(*it))
            return false;
    }

    return true;
}

/** Return whether or not any of the groups in this set contain any
    of the atoms of the view of the molecule in 'molview' */
bool MolGroupsBase::intersects(const MoleculeView &molview) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator
                                     it = molnum_to_mgnum.find(molview.data().number());

    if (it == molnum_to_mgnum.end())
        return false;

    QVarLengthArray<const MoleculeGroup*,10> groups;

    this->getGroups(*it, groups);

    int count = groups.count();
    const MoleculeGroup* const *groups_array = groups.constData();

    for (int i=0; i<count; ++i)
    {
        if (groups_array[i]->intersects(molview))
            return true;
    }

    return false;
}

/** Return whether any of the groups in this set contain any of the
    atoms of any of the views of any of the molecules in 'molecules' */
bool MolGroupsBase::intersects(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        if (this->intersects(*it))
            return true;
    }

    return false;
}

/** Return the list of molecule groups numbers of groups that
    contain at least one atom of the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
const QList<MGNum>& MolGroupsBase::groupsContaining(MolNum molnum) const
{
    QHash< MolNum,QList<MGNum> >::const_iterator it = molnum_to_mgnum.find(molnum);

    if (it == molnum_to_mgnum.end())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 is the groups in this set.")
                .arg(molnum), CODELOC );

    return *it;
}

/** Return the total number of groups in this set */
int MolGroupsBase::nGroups() const
{
    return mgidx_to_num.count();
}

/** Return the total number of groups in this set */
int MolGroupsBase::count() const
{
    return this->nGroups();
}

/** Return the total number of molecules in the groups in this set */
int MolGroupsBase::nMolecules() const
{
    return molnum_to_mgnum.count();
}

/** Return the total number of atoms in this groups in this set */
int MolGroupsBase::nAtoms() const
{
    int n = 0;

    for (const auto &mol : this->molecules())
    {
        n += mol.nAtoms();
    }

    return n;
}

/** Return the total number of residues in this groups in this set */
int MolGroupsBase::nResidues() const
{
    int n = 0;

    for (const auto &mol : this->molecules())
    {
        n += mol.nResidues();
    }

    return n;
}

/** Return the total number of chains in this groups in this set */
int MolGroupsBase::nChains() const
{
    int n = 0;

    for (const auto &mol : this->molecules())
    {
        n += mol.nChains();
    }

    return n;
}

/** Return the total number of segments in this groups in this set */
int MolGroupsBase::nSegments() const
{
    int n = 0;

    for (const auto &mol : this->molecules())
    {
        n += mol.nSegments();
    }

    return n;
}

/** Return the total number of views of molecules in the groups in this set.
    Note that if a view appears multiple times, then it will be counted
    multiple times */
int MolGroupsBase::nViews() const
{
    const QHash<MGNum,const MoleculeGroup*> groups = this->getGroups();

    int nviews = 0;

    for (QHash<MGNum,const MoleculeGroup*>::const_iterator it = groups.constBegin();
         it != groups.constEnd();
         ++it)
    {
        nviews += (*it)->nViews();
    }

    return nviews;
}

/** Return the total number of views of the molecule with number
    'molnum' in the groups in this set. If a view appears multiple
    times then it will be counted multiple times.

    \throw SireMol::missing_molecule
*/
int MolGroupsBase::nViews(MolNum molnum) const
{
    QVarLengthArray<const MoleculeGroup*,10> groups;
    this->getGroups(groupsContaining(molnum), groups);

    int count = groups.count();
    const MoleculeGroup* const *groups_array = groups.constData();

    int nviews = 0;

    for (int i=0; i<count; ++i)
    {
        nviews += groups_array[i]->nViews(molnum);
    }

    return nviews;
}

/** Return whether or not this set is empty (contains no groups) */
bool MolGroupsBase::isEmpty() const
{
    return mgidx_to_num.isEmpty();
}

/** Return the complete set of all molecules in this group. If a view of a
    molecule appears multiple times in this set then multiple copies of
    that view will be placed into the returned Molecules object.
    Note that this is a potentially very slow operation! */
Molecules MolGroupsBase::molecules() const
{
    Molecules all_mols;

    const QHash<MGNum,const MoleculeGroup*> groups = this->getGroups();

    for (QHash<MGNum,const MoleculeGroup*>::const_iterator it = groups.constBegin();
         it != groups.constEnd();
         ++it)
    {
        all_mols += (*it)->molecules();
    }

    return all_mols;
}

/** Return the complete set of all molecules in the group(s) that
    match the ID 'mgid'. If a view of a molecule appears multiple times
    in this set then multiple copies of that view will be placed into the
    returned molecules object.

    Note that this is potentially a very slow function
*/
Molecules MolGroupsBase::molecules(const MGID &mgid) const
{
    QList<MGNum> mgnums = mgid.map(*this);

    if (mgnums.count() == 1)
    {
        return this->at(mgnums.at(0)).molecules();
    }
    else
    {
        Molecules all_mols;

        foreach (MGNum mgnum, mgnums)
        {
            all_mols += this->at(mgnum).molecules();
        }

        return all_mols;
    }
}

/** Return the numbers of all molecule groups in this set */
QList<MGNum> MolGroupsBase::mgNums() const
{
    return mgidx_to_num;
}

/** Return the names of all molecule groups in this set */
QList<MGName> MolGroupsBase::mgNames() const
{
    QList<MGName> names;

    foreach (MGNum mgnum, mgidx_to_num)
    {
        names.append( this->operator[](mgnum).name() );
    }

    return names;
}

/** Assert that this set contains at least one atom of the
    molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
void MolGroupsBase::assertContains(MolNum molnum) const
{
    if (not molnum_to_mgnum.contains(molnum))
        throw SireMol::missing_molecule( QObject::tr(
            "None of the groups in this set contain the molecule with "
            "number %1.")
                .arg(molnum), CODELOC );
}

/** Assert that this set contains at least one atom of any
    molecule that is identified by the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
void MolGroupsBase::assertContains(const MolID &molid) const
{
    this->map(molid);
}

/** Assert that this contains the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
void MolGroupsBase::assertContains(MGNum mgnum) const
{
    if (not mgidx_to_num.contains(mgnum))
        throw SireMol::missing_group( QObject::tr(
            "This set does not contain the molecule group with "
            "number %1. Contained groups have numbers %2.")
                .arg(mgnum).arg(Sire::toString(mgidx_to_num)), CODELOC );
}

/** Assert that this contains at least one molecule group that
    is identified by the ID 'mgid'

    \throw SireMol:missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::assertContains(const MGID &mgid) const
{
    this->map(mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const MoleculeView &molview, const MGID &mgid)
{
    this->addIfUnique(molview, mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->addIfUnique(molviews, mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const Molecules &molecules, const MGID &mgid)
{
    this->addIfUnique(molecules, mgid);
}

/** Synonym for MolGroupsBase::addIfUnique

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MolGroupsBase::unite(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->addIfUnique(molgroup, mgid);
}

/** Remove the view of the molecule in 'molview' from all of
    the groups in this set. This does nothing if this exact
    view is not contained by any of the groups. If a group has
    multiple copies of this view, then this removes only
    the first copy. */
bool MolGroupsBase::remove(const MoleculeView &molview)
{
    //get the groups containing this molecule
    QList<MGNum> mgnums = molnum_to_mgnum.value( molview.data().number() );

    if (not mgnums.isEmpty())
        return this->remove( molview, IDOrSet<MGID>(mgnums) );
    else
        return false;
}

/** Remove the views of the molecule in 'molviews' from all of
    the groups in this set. This does nothing if none of these
    views are contained by any of the groups. If a group contains
    multiple copies of a view, then only the first copy is removed. */
bool MolGroupsBase::remove(const ViewsOfMol &molviews)
{
    //get the groups containing this molecules
    QList<MGNum> mgnums = molnum_to_mgnum.value( molviews.number() );

    if (not mgnums.isEmpty())
        return this->remove( molviews, IDOrSet<MGID>(mgnums) );
    else
        return false;
}

/** Remove all of the views of all of the molecules in 'molecules'
    from all of the groups in this set. If a group contains multiple
    copies of a view then only the first copy is removed */
bool MolGroupsBase::remove(const Molecules &molecules)
{
    if (molecules.isEmpty() or this->isEmpty())
        return false;

    return this->remove( molecules, IDOrSet<MGID>(mgidx_to_num) );
}

/** Remove all of the views of all of the molecules in 'molgroup'
    from this set. Note that this does not remove the actual molecule
    group. If you want to remove the group, then use the
    MolGroupsBase::remove(MGNum) function. Note that this also
    only removes the first copy of any duplicated views. */
bool MolGroupsBase::remove(const MoleculeGroup &molgroup)
{
    return this->remove(molgroup.molecules());
}

/** Remove all copies of the view of the molecule in 'molview'
    from all of the groups in this set. This removes all copies
    of a view (even duplicate copies) */
bool MolGroupsBase::removeAll(const MoleculeView &molview)
{
    QList<MGNum> mgnums = molnum_to_mgnum.value(molview.data().number());

    if (not mgnums.isEmpty())
        return this->removeAll( molview, IDOrSet<MGID>(mgnums) );
    else
        return false;
}

/** Remove all copies of the views of the molecule in 'molviews'
    from all of the groups in this set. This removes all
    copies of the views (even duplicate copies) */
bool MolGroupsBase::removeAll(const ViewsOfMol &molviews)
{
    QList<MGNum> mgnums = molnum_to_mgnum.value(molviews.number());

    if (not mgnums.isEmpty())
        return this->removeAll( molviews, IDOrSet<MGID>(mgnums) );
    else
        return false;
}

/** Remove all copies of all views of all molecules in 'molecules'
    from all of the groups in this set. This removes all copies
    of the views (even duplicate copies) */
bool MolGroupsBase::removeAll(const Molecules &molecules)
{
    if (molecules.isEmpty() or this->isEmpty())
        return false;

    return this->removeAll(molecules, IDOrSet<MGID>(mgidx_to_num));
}

/** Remove all copies of all views of all molecules in the
    group 'molgroup' from this set. Note that this removes
    the molecules, not the group. Note also that all copies
    of the views are removed (even duplicate copies) */
bool MolGroupsBase::removeAll(const MoleculeGroup &molgroup)
{
    return this->removeAll(molgroup.molecules());
}

/** Completely remove all views of the molecule with number
    'molnum' from all of the groups from this set. This
    does nothing if there are no views of this molecule
    in any of the groups  */
bool MolGroupsBase::remove(MolNum molnum)
{
    QList<MGNum> mgnums = molnum_to_mgnum.value(molnum);

    if (not mgnums.isEmpty())
        return this->remove(molnum, IDOrSet<MGID>(mgnums));
    else
        return false;
}

/** Completely remove all views of the molecules whose numbers
    are in 'molnums' from all of the groups in this set. This
    does nothing if there are no views of these molecules in
    any of the groups */
bool MolGroupsBase::remove(const QSet<MolNum> &molnums)
{
    if (molnums.isEmpty())
        return false;
    else if (molnums.count() == 1)
    {
        return this->remove( *(molnums.begin()) );
    }

    //get the list of groups that contain these molecules
    QList<MGNum> mgnums;

    foreach (MolNum molnum, molnums)
    {
        mgnums += molnum_to_mgnum.value(molnum);
    }

    mgnums = convert_to_qset(mgnums).values();

    if (not mgnums.isEmpty())
        return this->remove(molnums, IDOrSet<MGID>(mgnums));
    else
        return false;
}

/** Remove all molecules that match 'molid' from all groups */
bool MolGroupsBase::remove(const MolID &molid)
{
    try
    {
        const QList<MolNum> molnums = this->map(molid);

        return MolGroupsBase::remove( convert_to_qset(molnums) );
    }
    catch(...)
    {}

    return false;
}

/** Remove all groups (and molecules) that match the ID 'mgid' */
bool MolGroupsBase::remove(const MGID &mgid)
{
    return this->removeAll(mgid);
}

/** Completely clear all of the groups in this set */
bool MolGroupsBase::removeAll()
{
    if (this->nMolecules() > 1)
    {
        return this->removeAll( IDOrSet<MGID>(mgidx_to_num) );
    }
    else
        return false;
}

/** Update the copies in this set of the molecule viewed in 'molview'
    to use the same version as 'molview' */
void MolGroupsBase::update(const MoleculeView &molview, bool auto_commit)
{
    this->update(molview.data(), auto_commit);
}

/** Return a reference to the molecule data for the molecule whose data
    is in 'moldata' that is at the same version as the data that is
    already present in the groups of this set. This only returns
    'moldata' if this molecule is not in this set */
const MoleculeData&
MolGroupsBase::matchToExistingVersion(const MoleculeData &moldata) const
{
    MolNum molnum = moldata.number();

    QHash< MolNum,QList<MGNum> >::const_iterator
                                    it = molnum_to_mgnum.find(molnum);

    if (it == molnum_to_mgnum.end())
        //the molecule is not in this set - just return moldata
        return moldata;

    //get this molecule...
    const MoleculeData &existing_data = this->getGroup(it->first())
                                                 .molecule(molnum).data();

    return existing_data;
}

/** Return whether or not we need to update the molecule whose
    data is in 'moldata' */
bool MolGroupsBase::needToUpdate(const MoleculeData &moldata) const
{
    return moldata != this->matchToExistingVersion(moldata);
}

/** Return the set of molecules that is a copy of 'molecules', but where
    each molecule has been updated to match the version that exists in
    the groups in this set */
Molecules MolGroupsBase::matchToExistingVersion(const Molecules &molecules) const
{
    Molecules updated_mols = molecules;

    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        updated_mols.update( this->matchToExistingVersion(it->data()) );
    }

    return updated_mols;
}

/** Update the index with the new name of a molecule group */
void MolGroupsBase::changeNameIndex(MGNum mgnum, const MGName &old_name,
                                    const MGName &new_name)
{
    mgname_to_mgnum[old_name].removeAll(mgnum);

    if (mgname_to_mgnum[old_name].isEmpty())
        mgname_to_mgnum.remove(old_name);

    mgname_to_mgnum[new_name].append(mgnum);
}

/** Add the molecule group and all of its contained molecules to the index */
void MolGroupsBase::addToIndex(const MoleculeGroup &molgroup)
{
    if (mgidx_to_num.contains(molgroup.number()))
        //this group is already in the index!
        return;

    //record the number of the group - this allows the groups
    //to be indexed
    mgidx_to_num.append(molgroup.number());

    //record the name of the group - this allows the groups
    //to be accessed by name
    mgname_to_mgnum[molgroup.name()].append(molgroup.number());

    //now add the contents of this group to the index
    this->addToIndex(molgroup.number(), convert_to_qset(molgroup.molNums()));
}

/** Add the molecule with number 'molnum' to the index of the group
    with number 'mgnum' */
void MolGroupsBase::addToIndex(MGNum mgnum, MolNum molnum)
{
    BOOST_ASSERT( mgidx_to_num.contains(mgnum) );

    QHash< MolNum,QList<MGNum> >::const_iterator
                                it = molnum_to_mgnum.constFind(molnum);

    if (it != molnum_to_mgnum.constEnd() and
        it->contains(mgnum))
    {
        //this molecule is already recorded to be in the group
        return;
    }

    molnum_to_mgnum[molnum].append(mgnum);
}

/** Add the set of molecules whose numbers are in 'molnums' to the
    index of the group whose number is in 'mgnum' */
void MolGroupsBase::addToIndex(MGNum mgnum, const QSet<MolNum> &molnums)
{
    foreach (MolNum molnum, molnums)
    {
        this->addToIndex(mgnum, molnum);
    }
}

/** Add the set of molecules whose numbers are in 'molnums' to the
    index of the group whose number is in 'mgnum' */
void MolGroupsBase::addToIndex(MGNum mgnum, const QList<MolNum> &molnums)
{
    foreach (MolNum molnum, molnums)
    {
        this->addToIndex(mgnum, molnum);
    }
}

/** Completely remove the group with number 'mgnum' from the index */
void MolGroupsBase::removeFromIndex(MGNum mgnum)
{
    if (mgidx_to_num.removeAll(mgnum) == 0)
        //no groups were removed
        return;

    //ok, the group was removed - remove it from the other indicies
    this->clearIndex(mgnum);

    QMutableHashIterator< MGName,QList<MGNum> > it( mgname_to_mgnum );

    while (it.hasNext())
    {
        it.next();

        it.value().removeAll(mgnum);

        if (it.value().isEmpty())
            it.remove();
    }
}

/** Completely remove the molecule with number 'molnum' from the index */
void MolGroupsBase::removeFromIndex(MolNum molnum)
{
    molnum_to_mgnum.remove(molnum);
}

/** Remove the molecule with number 'molnum' from the index
    of the group with number 'mgnum' */
void MolGroupsBase::removeFromIndex(MGNum mgnum, MolNum molnum)
{
    QHash< MolNum,QList<MGNum> >::iterator it = molnum_to_mgnum.find(molnum);

    if (it != molnum_to_mgnum.end())
    {
        it->removeAll(mgnum);

        if (it->isEmpty())
            //there are no more groups that contain this molecule
            molnum_to_mgnum.remove(molnum);
    }
}

/** Remove all of the molecules whose numbers are in 'molnums' from the
    group with number 'mgnum' */
void MolGroupsBase::removeFromIndex(MGNum mgnum, const QSet<MolNum> &molnums)
{
    foreach (MolNum molnum, molnums)
    {
        this->removeFromIndex(mgnum, molnum);
    }
}

/** Completely clear the index of the group with number 'mgnum' - this
    removes the link between all molecules in this index and the group */
void MolGroupsBase::clearIndex(MGNum mgnum)
{
    QMutableHashIterator< MolNum,QList<MGNum> > it( molnum_to_mgnum );

    while (it.hasNext())
    {
        it.next();

        it.value().removeAll(mgnum);

        if (it.value().isEmpty())
            it.remove();
    }
}

/** Completely clear the entire index of molecules and molecule groups */
void MolGroupsBase::clearIndex()
{
    mgidx_to_num.clear();
    mgname_to_mgnum.clear();
    molnum_to_mgnum.clear();
}

int MolGroupsBase::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int MolGroupsBase::nFrames(const SireBase::PropertyMap &map) const
{
    return this->molecules().nFrames(map);
}

void MolGroupsBase::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void MolGroupsBase::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void MolGroupsBase::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void MolGroupsBase::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void MolGroupsBase::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    auto mols = this->molecules();
    mols.loadFrame(frame, map);
    this->update(mols);
}

void MolGroupsBase::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    auto mols = this->molecules();
    mols.saveFrame(frame, map);
    this->update(mols);
}

void MolGroupsBase::saveFrame(const SireBase::PropertyMap &map)
{
    auto mols = this->molecules();
    mols.saveFrame(map);
    this->update(mols);
}

void MolGroupsBase::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    auto mols = this->molecules();
    mols.deleteFrame(frame, map);
    this->update(mols);
}

/////////////
///////////// Implementation of MoleculeGroups
/////////////

RegisterMetaType<MoleculeGroups> r_molgroups;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const MoleculeGroups &molgroups)
{
    writeHeader(ds, r_molgroups, 1);

    SharedDataStream sds(ds);

    sds << molgroups.mgroups
        << static_cast<const MolGroupsBase&>(molgroups);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MoleculeGroups &molgroups)
{
    VersionID v = readHeader(ds, r_molgroups);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> molgroups.mgroups
            >> static_cast<MolGroupsBase&>(molgroups);
    }
    else
        throw version_error(v, "1", r_molgroups, CODELOC);

    return ds;
}

/** Construct an empty set of groups */
MoleculeGroups::MoleculeGroups()
          : ConcreteProperty<MoleculeGroups,MolGroupsBase>()
{}

const MoleculeGroups& MolGroupsBase::null()
{
    return *(create_shared_null<MoleculeGroups>());
}

/** Construct a set of groups that contains only the single group
    'molgroup' */
MoleculeGroups::MoleculeGroups(const MoleculeGroup &molgroup)
          : ConcreteProperty<MoleculeGroups,MolGroupsBase>()
{
    this->add(molgroup);
}

/** Construct a set of groups that contains the groups in 'molgroups'.
    The groups will be added in the order that they appear in this list. */
MoleculeGroups::MoleculeGroups(const QList<MolGroupPtr> &molgroups)
          : ConcreteProperty<MoleculeGroups,MolGroupsBase>()
{
    foreach (const MoleculeGroup &molgroup, molgroups)
    {
        this->add(molgroup);
    }
}

/** Copy constructor */
MoleculeGroups::MoleculeGroups(const MoleculeGroups &other)
          : ConcreteProperty<MoleculeGroups,MolGroupsBase>(other),
            mgroups(other.mgroups)
{}

/** Destructor */
MoleculeGroups::~MoleculeGroups()
{}

/** Copy assignment operator */
MoleculeGroups& MoleculeGroups::operator=(const MoleculeGroups &other)
{
    if (this != &other)
    {
        mgroups = other.mgroups;
        MolGroupsBase::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool MoleculeGroups::operator==(const MoleculeGroups &other) const
{
    return this == &other or mgroups == other.mgroups;
}

/** Comparison operator */
bool MoleculeGroups::operator!=(const MoleculeGroups &other) const
{
    return this != &other and mgroups != other.mgroups;
}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& MoleculeGroups::at(MGNum mgnum) const
{
    QHash<MGNum,MolGroupPtr>::const_iterator it = mgroups.find(mgnum);

    if (it == mgroups.end())
        throw SireMol::missing_group( QObject::tr(
            "There is no molecule group with number %1 in this set. "
            "The available groups are [ %2 ].")
                .arg(mgnum)
                .arg( Sire::toString(this->mgNums()) ), CODELOC );

    return it->read();
}

/** Completely reindex all of the groups */
void MoleculeGroups::reindex()
{
    //get the current order of groups
    QList<MGNum> mgnums = this->mgNums();

    //completely clear the index
    MolGroupsBase::clearIndex();

    auto remaining_groups = convert_to_qset(mgroups.keys());

    //reindex the groups, in order
    foreach (MGNum mgnum, mgnums)
    {
        if (remaining_groups.contains(mgnum))
        {
            this->addToIndex( mgroups.value(mgnum) );
            remaining_groups.remove(mgnum);
        }
    }

    //add any remaining groups
    foreach (MGNum mgnum, remaining_groups)
    {
        this->addToIndex( mgroups.value(mgnum) );
    }
}

/** Add the molecule group 'molgroup' to this set. This does
    nothing if this group is already in this set. This updates
    the molecules in 'molgroup' so that they are at the
    same version as any existing copies of the molecules
    in this set. */
void MoleculeGroups::add(const MoleculeGroup &molgroup)
{
    if (this->contains(molgroup.number()))
        return;

    //copy this group
    MolGroupPtr new_group(molgroup);

    //update the group so that it has the same version of the
    //molecules as the other groups in this set
    for (QHash<MGNum,MolGroupPtr>::const_iterator it = mgroups.constBegin();
         it != mgroups.constEnd();
         ++it)
    {
        new_group.edit().update(it->read().molecules());
    }

    //add the group to the index
    this->addToIndex(new_group.read());

    mgroups.insert(new_group.read().number(), new_group);
}

/** Addition operator */
MoleculeGroups& MoleculeGroups::operator+=(const MoleculeGroup &molgroup)
{
    this->add(molgroup);
    return *this;
}

/** Update the group 'molgroup'. If this group is in this set,
    then it updates the group to the same version. Then, regardless
    of whether the group is in this set, it then updates all
    molecules in all of the groups so that they have the same
    version number as 'molgroup'. This does nothing if
    molgroup and none of its molecules are in this set */
void MoleculeGroups::update(const MoleculeGroup &molgroup, bool auto_commit)
{
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        copy.accept();
        this->update(copy);
        return;
    }

    //do this in a copy, as something weird may go wrong...
    MoleculeGroups orig_groups(*this);

    if (this->needsAccepting())
        this->accept();

    try
    {

    QHash<MGNum,MolGroupPtr>::iterator it = mgroups.find(molgroup.number());

    if (it != mgroups.end() and
        molgroup.majorVersion() != it->read().majorVersion())
    {
        //this group exists in this set and it is at a different
        //major version - this means that its complement of
        //molecules has changed, so there needs to be a change
        //in the index...

        //get the list of current molecules in this group
        auto old_molnums = convert_to_qset(it->read().molNums());

        //now the new set of numbers...
        auto new_molnums = convert_to_qset(molgroup.molNums());

        this->addToIndex(molgroup.number(), new_molnums - old_molnums);
        this->removeFromIndex(molgroup.number(), old_molnums - new_molnums);
    }

    //update all of the groups;
    for (it = mgroups.begin() ; it != mgroups.end(); ++it)
    {
        it->edit().update(molgroup);
    }

    }
    catch(...)
    {
        //something went wrong - restore the original...
        this->operator=(orig_groups);
        throw;
    }
}

/** Completely remove the group with number 'mgnum'. This
    does nothing if there is no such group in this set */
bool MoleculeGroups::_pvt_remove(MGNum mgnum)
{
    if (mgroups.contains(mgnum))
    {
        mgroups.remove(mgnum);
        this->removeFromIndex(mgnum);
        return true;
    }

    return false;
}

/** Remove the molecules contained in 'molgroup' from this set.
    Note that this *does not* remove this molecule group itself
     - if you want to remove the molecule group, use
     MoleculeGroups::remove(molgroup.number()) */
bool MoleculeGroups::remove(const MoleculeGroup &molgroup)
{
    return MolGroupsBase::remove(molgroup.molecules());
}

/** Remove the groups that match the ID 'mgid' from this set. This
    does nothing if there are no such groups. */
bool MoleculeGroups::remove(const MGID &mgid)
{
    try
    {
        QList<MGNum> mgnums = this->map(mgid);

        bool removed = false;

        foreach (MGNum mgnum, mgnums)
        {
            if (this->_pvt_remove(mgnum))
            {
                removed = true;
            }
        }

        return removed;
    }
    catch(...)
    {}

    return false;
}

/** Remove the molecules that match the ID 'molid' from this set.
    This does nothing if there are no molecules that match this
    ID in this set */
bool MoleculeGroups::remove(const MolID &molid)
{
    try
    {
        return MolGroupsBase::remove(convert_to_qset(this->map(molid)));
    }
    catch(...)
    {}

    return false;
}

/** Remove operator */
MoleculeGroups& MoleculeGroups::operator-=(const MolID &molid)
{
    this->remove(molid);
    return *this;
}

/** Remove operator */
MoleculeGroups& MoleculeGroups::operator-=(const MGID &mgid)
{
    this->remove(mgid);
    return *this;
}

/** Remove operator */
MoleculeGroups& MoleculeGroups::operator-=(const MoleculeGroup &molgroup)
{
    this->remove(molgroup);
    return *this;
}

/** Remove operator */
MoleculeGroups& MoleculeGroups::operator-=(const Molecules &molecules)
{
    this->remove(molecules);
    return *this;
}

/** Add the view of the molecule in 'molview' to the groups
    identified by 'mgid'. This adds the view as a duplicate
    if it already exists in the group. The version
    of the molecule added is the version already present
    in this set, if it exists.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::add(const MoleculeView &molview, const MGID &mgid)
{
    //get the number of the groups identified by 'mgid'
    QList<MGNum> mgnums = this->map(mgid);

    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->edit().add(view);
        this->addToIndex(mgnum, view.number());
    }
}

/** Add the views of the molecule in 'molviews' to the groups
    identified by 'mgid'. This adds the view as a duplicate if
    it already exists in a group.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::add(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->edit().add(views);
        this->addToIndex(mgnum, views.number());
    }
}

/** Add each of the molecules in 'molecules' to the groups
    identified by 'mgid'. This adds the views as duplicates
    if they exist already in a group. Any molecules that
    already exist in any of the groups in this set are
    updated to the versions that are already present
    in this set.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::add(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    Molecules mols = this->matchToExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();

    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->edit().add(mols);
        this->addToIndex(mgnum, molnums);
    }
}

/** Add the molecules in the group 'molgroup' to the groups
    identified by 'mgid'. This adds the views as duplicates
    if they already exist, and adds the views in the same
    order as they appear in 'molgroup'. This is slightly less
    efficient than MoleculeGroups::add(const Molecules&), so use
    that function if you don't care about the order.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::add(const MoleculeGroup &molgroup, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    //update the group to match the molecule versions
    //already present in this group...
    MolGroupPtr group(molgroup);
    group.edit().update( this->matchToExistingVersion(group.read().molecules()) );

    auto molnums = convert_to_qset(group.read().molNums());

    foreach (MGNum mgnum, mgnums)
    {
        mgroups.find(mgnum)->edit().add(group);
        this->addToIndex(mgnum,molnums);
    }
}

/** Add the view of the molecule in 'molview' to the groups
    identified by 'mgid'. This only adds the view to a group
    if it doesn't already exist in the group. The version
    of the molecule already present in this set is used if
    such a molecule already exists.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::addIfUnique(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        if ( mgroups.find(mgnum)->edit().addIfUnique(view) )
            this->addToIndex(mgnum,view.number());
    }
}

/** Add the views of the molecule in 'molviews' ot the groups
    identified by 'mgid'. This only adds views to groups that
    don't already exist in that group, and uses the existing
    version of the molecule is it already exists in one
    of the groups of this set.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        ViewsOfMol added_views = mgroups.find(mgnum)->edit().addIfUnique(views);

        if (not added_views.isEmpty())
            this->addToIndex(mgnum,views.number());
    }
}

/** Add all of the views of the molecules in 'molecules' to the groups
    identified by 'mgid'. This only adds views that don't already
    exist in the group, and uses the version of the molecules that already
    exists in one of the groups of this set (if one exists)

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::addIfUnique(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    Molecules mols = this->matchToExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();

    foreach (MGNum mgnum, mgnums)
    {
        QList<ViewsOfMol> added_mols = mgroups.find(mgnum)->edit().addIfUnique(mols);

        if (not added_mols.isEmpty())
        {
            QSet<MolNum> added_molnums;

            foreach (const ViewsOfMol &mol, added_mols)
            {
                added_molnums.insert(mol.number());
            }

            this->addToIndex(mgnum,added_molnums);
        }
    }
}

/** This adds all of the views of the molecules in the group
    'molgroup', in the same order as they exist in this group,
    to all of the groups identified by 'mgid'. This only
    adds views to a group that don't already exist in that
    group and uses the existing version of the molecule if
    it exists anywhere in this set.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    //update the group...
    MoleculeGroup group(molgroup);
    group.update( this->matchToExistingVersion(group.molecules()) );

    //now add the molecules...
    foreach (MGNum mgnum, mgnums)
    {
        QList<ViewsOfMol> added_mols = mgroups.find(mgnum)->edit().addIfUnique(group);

        if (not added_mols.isEmpty())
        {
            QSet<MolNum> added_molnums;

            foreach (const ViewsOfMol &mol, added_mols)
            {
                added_molnums.insert(mol.number());
            }

            this->addToIndex(mgnum,added_molnums);
        }
    }
}

/** Remove the view of the molecule in 'molview' from the groups
    identified by 'mgid'. This only removes the first copy
    of the view from each group, if multiple copies exist.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::remove(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    MolNum molnum = molview.data().number();

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        if (mgroup.edit().remove(molview))
        {
            removed_mol = true;

            if (not mgroup.read().contains(molnum))
                this->removeFromIndex(mgnum,molnum);
        }
    }

    return removed_mol;
}

/** Remove the views of the molecule in 'molviews' from the groups
    identified by 'mgid'. This only removes the first copy of the
    views from each group if they exist multiple times.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::remove(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    MolNum molnum = molviews.number();

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        ViewsOfMol removed_views = mgroup.edit().remove(molviews);

        if (not removed_views.isEmpty())
        {
            removed_mol = true;

            if (not mgroup.read().contains(molnum))
                this->removeFromIndex(mgnum,molnum);
        }
    }

    return removed_mol;
}

/** Remove all of the views of the molecules in 'molecules' from
    the groups identified by 'mgid'. This removes only the first
    copies of the views in each group.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::remove(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        QList<ViewsOfMol> removed_mols = mgroup.edit().remove(molecules);

        if (not removed_mols.isEmpty())
            removed_mol = true;

        QSet<MolNum> removed_molnums;

        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            if (not mgroup.read().contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }

        if (not removed_molnums.isEmpty())
            this->removeFromIndex(mgnum, removed_molnums);
    }

    return removed_mol;
}

/** Remove all of the views of the molecules in the group 'molgroup' from
    the groups identified by 'mgid'. This removes only the first
    copies of the views in each group.

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::remove(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->remove(molgroup.molecules(), mgid);
}

/** Remove all copies of the view of the molecule in 'molview' from
    the groups identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::removeAll(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    MolNum molnum = molview.data().number();

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        if (mgroup.edit().removeAll(molview))
        {
            removed_mol = true;

            if (not mgroup.read().contains(molnum))
                this->removeFromIndex(mgnum, molnum);
        }
    }

    return removed_mol;
}

/** Remove all copies of the views of the molecule in 'molviews' from
    the groups identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::removeAll(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    MolNum molnum = molviews.number();

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        ViewsOfMol removed_views = mgroup.edit().removeAll(molviews);

        if (not removed_views.isEmpty())
        {
            removed_mol = true;

            if (not mgroup.read().contains(molnum))
                this->removeFromIndex(mgnum, molnum);
        }
    }

    return removed_mol;
}

/** Remove all copies of the views of the molecules in 'molecules'
    from the groups identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::removeAll(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        QList<ViewsOfMol> removed_mols = mgroup.edit().removeAll(molecules);

        if (not removed_mols.isEmpty())
            removed_mol = true;

        QSet<MolNum> removed_molnums;

        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            if (not mgroup.read().contains(removed_mol.number()))
                removed_molnums.insert(removed_mol.number());
        }

        this->removeFromIndex(mgnum,removed_molnums);
    }

    return removed_mol;
}

/** Remove all of the molecules from all of the groups identified by
    the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::removeAll(const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MoleculeGroup &molgroup = mgroups.find(mgnum)->edit();

        if (not molgroup.isEmpty())
        {
            removed_mol = true;
            molgroup.removeAll();
            this->removeFromIndex(mgnum);
        }
    }

    return removed_mol;
}

/** Remove all of the views of the molecules in the group 'molgroup'
    from the groups identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::removeAll(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->removeAll(molgroup.molecules(), mgid);
}

/** Remove all views of the molecule with number 'molnum' from the
    groups identified by 'mgid'. This does nothing to any groups
    that don't contain this molecule

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::remove(MolNum molnum, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        ViewsOfMol removed_views = mgroup.edit().remove(molnum);

        if (not removed_views.isEmpty())
        {
            removed_mol = true;
            this->removeFromIndex(mgnum, molnum);
        }
    }

    return removed_mol;
}

/** Remove the molecules whose numbers are in 'molnums' from the
    groups identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool MoleculeGroups::remove(const QSet<MolNum> &molnums, const MGID &mgid)
{
    QList<MGNum> mgnums = this->map(mgid);

    bool removed_mol = false;

    foreach (MGNum mgnum, mgnums)
    {
        MolGroupPtr &mgroup = *(mgroups.find(mgnum));

        QList<ViewsOfMol> removed_mols = mgroup.edit().remove(molnums);

        if (not removed_mols.isEmpty())
            removed_mol = true;

        QSet<MolNum> removed_nums;

        foreach (const ViewsOfMol &removed_mol, removed_mols)
        {
            removed_nums.insert(removed_mol.number());
        }

        this->removeFromIndex(mgnum, removed_nums);
    }

    return removed_mol;
}

/** Update all of the groups to use the version of the molecule
    present in 'moldata' */
void MoleculeGroups::update(const MoleculeData &moldata, bool auto_commit)
{
    if (auto_commit and this->needsAccepting())
    {
        this->accept();
    }

    //get the current copy of the molecule...
    if (this->needToUpdate(moldata))
    {
        const QList<MGNum> &mgnums = this->groupsContaining(moldata.number());

        foreach (MGNum mgnum, mgnums)
        {
            mgroups.find(mgnum)->edit().update(moldata, auto_commit);
        }
    }

    if (auto_commit and this->needsAccepting())
    {
        this->accept();
    }
}

/** Update all of the groups to use the versions of the molecules
    held in 'molecules' */
void MoleculeGroups::update(const Molecules &molecules, bool auto_commit)
{
    //we need to do this in a copy, so that we can revert
    //if an error occurs...
    boost::shared_ptr<MoleculeGroups> old_state( this->clone() );

    try
    {
        if (auto_commit and this->needsAccepting())
            this->accept();

        for (QHash<MGNum,MolGroupPtr>::iterator it = mgroups.begin();
             it != mgroups.end();
             ++it)
        {
            it->edit().update(molecules, auto_commit);
        }

        if (auto_commit and this->needsAccepting())
            this->accept();
    }
    catch(...)
    {
        //something went wrong...
        this->operator=(*old_state);
        throw;
    }
}

/** Return whether or not this set of molecule groups is using a temporary
    workspace and needs accepting */
bool MoleculeGroups::needsAccepting() const
{
    for (QHash<MGNum,MolGroupPtr>::const_iterator it = mgroups.constBegin();
         it != mgroups.constEnd();
         ++it)
    {
        if (it->read().needsAccepting())
            return true;
    }

    return false;
}

/** Tell the molecule group that the last move was accepted. This tells the
    group to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void MoleculeGroups::accept()
{
    if (needsAccepting())
    {
        for (QHash<MGNum,MolGroupPtr>::iterator it = mgroups.begin();
             it != mgroups.end();
             ++it)
        {
            if (it->read().needsAccepting())
                it->edit().accept();
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they only contain the view in 'molview'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::setContents(const MGID &mgid, const MoleculeView &molview)
{
    QList<MGNum> mgnums = this->map(mgid);

    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        if (mgroups.find(mgnum)->edit().setContents(view))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, view.number());
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they only contain the views of the molecule in 'molviews'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::setContents(const MGID &mgid, const ViewsOfMol &molviews)
{
    QList<MGNum> mgnums = this->map(mgid);

    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );

    foreach (MGNum mgnum, mgnums)
    {
        if (mgroups.find(mgnum)->edit().setContents(views))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, views.number());
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they contain only the views of the molecules contained in 'molecules'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::setContents(const MGID &mgid, const Molecules &molecules)
{
    QList<MGNum> mgnums = this->map(mgid);

    Molecules mols = this->matchToExistingVersion(molecules);

    QSet<MolNum> molnums = mols.molNums();

    foreach (MGNum mgnum, mgnums)
    {
        if (mgroups.find(mgnum)->edit().setContents(mols))
        {
            this->clearIndex(mgnum);
            this->addToIndex(mgnum,molnums);
        }
    }
}

/** Set the contents of the groups identified by 'mgid' so that
    they contain the same views of the same molecules in the
    same order as in the group 'molgroup'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
void MoleculeGroups::setContents(const MGID &mgid, const MoleculeGroup &molgroup)
{
    QList<MGNum> mgnums = this->map(mgid);

    QHash<MGNum,MolGroupPtr>::const_iterator
                                  it = mgroups.constFind(molgroup.number());

    if (it != mgroups.constEnd() and
        it->read().majorVersion() == molgroup.majorVersion() and
        it->read().minorVersion() == molgroup.minorVersion())
    {
        //the group is the same as the existing one, and at the same
        //version, so nothing to do :-)
        return;
    }

    //ok we need to set the contents, but not necessarily at the same
    //versions... (we need to update this group to use the versions
    //of any molecules already present in this group)
    Molecules mols = this->matchToExistingVersion(molgroup.molecules());

    QSet<MolNum> molnums = mols.molNums();

    foreach (MGNum mgnum, mgnums)
    {
        //try to set the contents using the wrong molecule version...
        MolGroupPtr &group = *(mgroups.find(mgnum));

        if (group.edit().setContents(molgroup))
        {
            //ok, we need to update the group to ensure that it uses
            //the same version of molecules as are present in the
            //rest of the set
            group.edit().update(mols);

            //now update the index
            this->clearIndex(mgnum);
            this->addToIndex(mgnum, molnums);
        }
    }
}

/** Protected function used to return a const reference to the
    group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& MoleculeGroups::getGroup(MGNum mgnum) const
{
    QHash<MGNum,MolGroupPtr>::const_iterator it = mgroups.find(mgnum);

    if (it == mgroups.end())
        throw SireMol::missing_group( QObject::tr(
            "Cannot find the MoleculeGroup with number %1. Available "
            "groups are %2")
                .arg(mgnum).arg(Sire::toString(mgroups.keys())), CODELOC );

    return it->read();
}

/** Protected function used to return const pointers to the
    groups whose numbers are in 'mgnums'

    \throw SireMol::missing_group
*/
void MoleculeGroups::getGroups(const QList<MGNum> &mgnums,
                          QVarLengthArray<const MoleculeGroup*,10> &groups) const
{
    groups.clear();

    QHash<MGNum,MolGroupPtr>::const_iterator it;

    foreach (MGNum mgnum, mgnums)
    {
        it = mgroups.find(mgnum);

        if (it == mgroups.end())
            throw SireMol::missing_group( QObject::tr(
                "Cannot find the MoleculeGroup with number %1. Available "
                "groups are %2")
                    .arg(mgnum).arg(Sire::toString(mgroups.keys())), CODELOC );


        groups.append( it->constData() );
    }
}

/** Protected function used to return a hash of const
    pointers to all of the groups in this set. */
QHash<MGNum,const MoleculeGroup*> MoleculeGroups::getGroups() const
{
    QHash<MGNum,const MoleculeGroup*> groups;

    groups.reserve(mgroups.count());

    for (QHash<MGNum,MolGroupPtr>::const_iterator it = mgroups.begin();
         it != mgroups.end();
         ++it)
    {
        groups.insert( it.key(), it->constData() );
    }

    return groups;
}

const char* MoleculeGroups::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MoleculeGroups>() );
}

MoleculeGroups* MoleculeGroups::clone() const
{
    return new MoleculeGroups(*this);
}
