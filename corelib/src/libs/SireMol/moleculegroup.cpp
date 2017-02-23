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

#include <boost/tuple/tuple.hpp>

#include "molnum.h"
#include "SireID/index.h"

namespace SireMol
{
static bool operator==(const boost::tuple<SireMol::MolNum,SireID::Index> &i0,
                       const boost::tuple<SireMol::MolNum,SireID::Index> &i1)
{
    return i0.get<0>() == i1.get<0>() and
           i0.get<1>() == i1.get<1>();
}
}

#include <QMutex>
#include <QVector>

#include "moleculegroup.h"
#include "partialmolecule.h"
#include "molecule.h"

#include "molnum.h"
#include "molname.h"
#include "molidx.h"
#include "molidentifier.h"

#include "mover.hpp"
#include "editor.hpp"

#include "mgname.h"
#include "mgnum.h"

#include "tostring.h"

#include "SireBase/incremint.h"
#include "SireBase/majorminorversion.h"
#include "SireBase/refcountdata.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using boost::tuple;

using namespace SireMol;
using namespace SireID;
using namespace SireStream;
using namespace SireBase;

////////
//////// Objects relating to the global registry of MoleculeGroups
////////

namespace SireBase
{
    template class VersionRegistry<MGNum>;
}

static VersionRegistry<MGNum> mgnum_registry;

static Incremint last_number(0);

MGNum MGNum::getUniqueNumber()
{
    MGNum new_num( last_number.increment() );

    while ( mgnum_registry.registered(new_num) )
    {
        new_num = MGNum(last_number.increment());
    }
    
    return new_num;
}

////////
//////// Functions relating to the tuple<MolNum,Index>
////////

static QDataStream& operator<<(QDataStream &ds, 
                               const tuple<MolNum,Index> &molviewidx)
{
    ds << molviewidx.get<0>() << molviewidx.get<1>();
    return ds;
}

static QDataStream& operator>>(QDataStream &ds, 
                               tuple<MolNum,Index> &molviewidx)
{
    MolNum molnum;
    Index viewidx;
    
    ds >> molnum >> viewidx;

    molviewidx = tuple<MolNum,Index>(molnum,viewidx);
    
    return ds;
}

////////////
//////////// Implementation of MolGroupPvt
////////////

namespace SireMol
{
namespace detail
{

class MolGroupPvt : public RefCountData
{
public:
    MolGroupPvt();
    
    MolGroupPvt(const QString &name);
    MolGroupPvt(const QString &name, const MolGroupPvt &other);
    
    MolGroupPvt(const MolGroupPvt &other);
    
    ~MolGroupPvt();

    MolGroupPvt& operator=(const MolGroupPvt &other);
    
    bool operator==(const MolGroupPvt &other) const;
    bool operator!=(const MolGroupPvt &other) const;
    
    void incrementMajor();
    void incrementMinor();
    
    Molecules molecules;
    
    QVector<MolNum> molidx_to_num;
    QVector< tuple<MolNum,Index> > molviewidx_to_num;
    
    MGName name;
    MGNum number;
    
    MajorMinorVersion version;
};

} // end of namespace detail
} // end of namespace SireMol;

using namespace SireMol::detail;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const MolGroupPvt &molgrouppvt)
{
    SharedDataStream sds(ds);
    
    QHash< MolNum,QList<MolNum> > molname_to_num;
    
    sds << molgrouppvt.molecules
        << molgrouppvt.molidx_to_num
        << molgrouppvt.molviewidx_to_num
        << molname_to_num
        << molgrouppvt.name
        << molgrouppvt.number;

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,     
                                       MolGroupPvt &molgrouppvt)
{
    SharedDataStream sds(ds);

    QHash< MolNum,QList<MolNum> > molname_to_num;
       
    sds >> molgrouppvt.molecules
        >> molgrouppvt.molidx_to_num
        >> molgrouppvt.molviewidx_to_num
        >> molname_to_num
        >> molgrouppvt.name
        >> molgrouppvt.number;

    molgrouppvt.version = mgnum_registry.registerObject(molgrouppvt.number);
    
    return ds;
}

static SharedDataPointer<MolGroupPvt> shared_null;

static const SharedDataPointer<MolGroupPvt>& getSharedNull()
{
    if (shared_null.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null.constData() == 0)
            shared_null = new MolGroupPvt();
    }
     
    return shared_null;
}

/** Construct an empty, unnamed MolGroupPvt */
MolGroupPvt::MolGroupPvt() 
            : RefCountData(),
              name("unnamed"),
              number( MGNum::getUniqueNumber() ),
              version(1,0)
{}

/** Construct an empty, named MolGroupPvt */
MolGroupPvt::MolGroupPvt(const QString &nme)
            : RefCountData(),
              name(nme),
              number( MGNum::getUniqueNumber() ),
              version(1,0)
{}

/** Construct a named group that contains the same molecules as 'other' */
MolGroupPvt::MolGroupPvt(const QString &nme, 
                         const MolGroupPvt &other)
            : RefCountData(),
              molecules(other.molecules),
              molidx_to_num(other.molidx_to_num),
              molviewidx_to_num(other.molviewidx_to_num),
              name(nme),
              number( MGNum::getUniqueNumber() ),
              version(1,0)
{}

/** Copy constructor */
MolGroupPvt::MolGroupPvt(const MolGroupPvt &other)
            : RefCountData(),
              molecules(other.molecules),
              molidx_to_num(other.molidx_to_num),
              molviewidx_to_num(other.molviewidx_to_num),
              name(other.name),
              number(other.number),
              version(other.version)
{}
                   
/** Destructor */
MolGroupPvt::~MolGroupPvt()
{}

/** Copy assignment operator */
MolGroupPvt& MolGroupPvt::operator=(const MolGroupPvt &other)
{
    if (this != &other)
    {
        molecules = other.molecules;
        molidx_to_num = other.molidx_to_num;
        molviewidx_to_num = other.molviewidx_to_num;
        name = other.name;
        number = other.number;
        version = other.version;
    }
    
    return *this;
}

bool MolGroupPvt::operator==(const MolGroupPvt &other) const
{
    return number == other.number and
           version == other.version;
}

bool MolGroupPvt::operator!=(const MolGroupPvt &other) const
{
    return number != other.number or
           version != other.version;
}

void MolGroupPvt::incrementMajor()
{
    version.incrementMajor();
}

void MolGroupPvt::incrementMinor()
{
    version.incrementMinor();
}

////////////
//////////// Implementation of MoleculeGroup
////////////

static const RegisterMetaType<MoleculeGroup> r_MoleculeGroup;

/** Serialise a MoleculeGroup to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const MoleculeGroup &molgroup)
{
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        copy.accept();
        ds << copy;
    }
    else
    {
        writeHeader(ds, r_MoleculeGroup, 1);
        SharedDataStream sds(ds);
        sds << molgroup.d;
    }
    
    return ds;
}

/** Deserialise a MoleculeGroup from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       MoleculeGroup &molgroup)
{
    VersionID v = readHeader(ds, r_MoleculeGroup);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        molgroup.workspace.clear();
        sds >> molgroup.d;
    }
    else
        throw version_error(v, "1", r_MoleculeGroup, CODELOC);
        
    return ds;
}

/** Default constructor */
MoleculeGroup::MoleculeGroup() 
              : ConcreteProperty<MoleculeGroup,Property>(),
                d( getSharedNull() )
{}


static SharedPolyPointer<MoleculeGroup> shared_null_molgroup;

const MoleculeGroup& MoleculeGroup::null()
{
    if (shared_null_molgroup.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null_molgroup.constData() == 0)
            shared_null_molgroup = new MoleculeGroup();
    }

    return *(shared_null_molgroup.constData());
}

/** Construct a group that holds the passed molecules */
MoleculeGroup::MoleculeGroup(const Molecules &molecules)
              : ConcreteProperty<MoleculeGroup,Property>(),
                d( getSharedNull() )
{
    this->add(molecules);
}

/** Construct an empty, but named, group */
MoleculeGroup::MoleculeGroup(const QString &name)
              : ConcreteProperty<MoleculeGroup,Property>(),
                d( new MolGroupPvt(name) )
{}

/** Construct a named group that contains the passed molecule */
MoleculeGroup::MoleculeGroup(const QString &name, const MoleculeView &molecule)
              : ConcreteProperty<MoleculeGroup,Property>(),
                d( new MolGroupPvt(name) )
{
    this->add(molecule);
}

/** Construct a named group that contains the passed molecules */
MoleculeGroup::MoleculeGroup(const QString &name, const Molecules &molecules)
              : ConcreteProperty<MoleculeGroup,Property>(),
                d( new MolGroupPvt(name) )
{
    this->add(molecules);
}
  
/** Construct a named group that contains the same molecules as 'other' */       
MoleculeGroup::MoleculeGroup(const QString &name, const MoleculeGroup &other)
              : ConcreteProperty<MoleculeGroup,Property>(), workspace(other.workspace)
{
    if (name == other.name())
    {
        d = other.d;
    }
    else
        d = new MolGroupPvt(name, *(other.d));
}

/** Copy constructor */
MoleculeGroup::MoleculeGroup(const MoleculeGroup &other)
              : ConcreteProperty<MoleculeGroup,Property>(),
                d(other.d), workspace(other.workspace)
{}

/** Destructor */
MoleculeGroup::~MoleculeGroup()
{}

/** Copy assignment operator */
MoleculeGroup& MoleculeGroup::operator=(const MoleculeGroup &other)
{
    Property::operator=(other);
    d = other.d;
    workspace = other.workspace;
    return *this;
}

/** Comparison operator */
bool MoleculeGroup::operator==(const MoleculeGroup &other) const
{
    return (d == other.d or *d == *(other.d)) and workspace == other.workspace;
}

/** Comparison operator */
bool MoleculeGroup::operator!=(const MoleculeGroup &other) const
{
    return not operator==(other);
}

/** Return the views of the molecule with number 'molnum' from this group.

    \throw SireMol::missing_molecule
*/
const ViewsOfMol& MoleculeGroup::at(MolNum molnum) const
{
    const ViewsOfMol &oldmol = d->molecules.at(molnum);
    
    if (workspace.isEmpty())
    {
        return oldmol;
    }
    else
        return workspace.getUpdated(oldmol);
}

/** Return the specified view of the specified molecule...

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
PartialMolecule MoleculeGroup::operator[](const boost::tuple<MolNum,Index> &viewidx) const
{
    PartialMolecule oldmol = d->molecules.at(viewidx);
    
    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the views of the molecule with number 'molnum' from this group

    \throw SireMol::missing_molecule
*/
const ViewsOfMol& MoleculeGroup::operator[](MolNum molnum) const
{
    return this->at(molnum);
}

/** Return the views of the molecule at index 'molidx'

    \throw SireMol::invalid_index
*/
const ViewsOfMol& MoleculeGroup::operator[](MolIdx molidx) const
{
    return this->at(molidx);
}

/** Return the views of the molecule called 'molname'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
const ViewsOfMol& MoleculeGroup::operator[](const MolName &molname) const
{
    return this->at(molname);
}

/** Return the views of the molecule that matches the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::operator[](const MolID &molid) const
{
    return this->at(molid);
}

/** Return the specified view of the specified molecule...

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
PartialMolecule 
MoleculeGroup::operator[](const boost::tuple<MolIdentifier,Index> &viewidx) const
{
    return this->at(viewidx);
}

/** Add the views of the molecules in 'molecules' to this group. This
    adds duplicates of any views that already exist in this group. */
MoleculeGroup& MoleculeGroup::operator+=(const Molecules &molecules)
{
    accept();
    this->add(molecules);
    return *this;
}

/** Remove the views of the molecules in 'molecules' from this
    group. This only removes the first copy of any duplicated
    views that exist in this group. */
MoleculeGroup& MoleculeGroup::operator-=(const Molecules &molecules)
{
    accept();
    this->remove(molecules);
    return *this;
}

/** Return a string representation of this MoleculeGroup */
QString MoleculeGroup::toString() const
{
    return QString("MoleculeGroup(%1, %2)").arg(this->name()).arg(this->number());
}

/** Return the version number of the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
quint64 MoleculeGroup::getMoleculeVersion(MolNum molnum) const
{
    return this->at(molnum).version();
}

/** Return the version number of the molecule with ID 'molid' 

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
quint64 MoleculeGroup::getMoleculeVersion(const MolID &molid) const
{
    return this->getMoleculeVersion( this->getMoleculeNumber(molid) );
}

/** Return the views of the molecule at index 'molidx'

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::at(MolIdx molidx) const
{
    return this->at( this->getMoleculeNumber(molidx) );
}

/** Return the views of the molecule called 'molname'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireMol::invalid_index
*/
const ViewsOfMol& MoleculeGroup::at(const MolName &molname) const
{
    return this->at( this->getMoleculeNumber(molname) );
}

/** Return the views of the molecule that is identified by 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::at(const MolID &molid) const
{
    return this->at( this->getMoleculeNumber(molid) );
}

/** Return the specified view of the specified molecule in this group.

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
PartialMolecule MoleculeGroup::at(const boost::tuple<MolNum,Index> &viewidx) const
{
    PartialMolecule oldmol = d->molecules.at(viewidx);
    
    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the view of hte molecule at 'viewidx'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
PartialMolecule 
MoleculeGroup::at(const boost::tuple<MolIdentifier,Index> &viewidx) const
{
    PartialMolecule oldmol = this->at( viewidx.get<0>().base(), viewidx.get<1>() );

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the specified view of the specified molecule in this group.

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
PartialMolecule MoleculeGroup::at(MolNum molnum, int viewidx) const
{
    PartialMolecule oldmol = d->molecules.at(molnum, viewidx);

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the specified view of the molecule identified by 
    the ID 'molid' 
    
    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
PartialMolecule MoleculeGroup::at(const MolID &molid, int viewidx) const
{
    PartialMolecule oldmol = d->molecules.at( this->getMoleculeNumber(molid),
                                              viewidx );

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the views of the molecule at index 'idx' in this group.

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::moleculeAt(int idx) const
{
    const ViewsOfMol &oldmol = d->molecules.at( d->molidx_to_num.constData()
                                        [Index(idx).map(d->molidx_to_num.count())] );

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the view of the molecule at index 'idx' in this group.

    \throw SireError::invalid_index
*/
PartialMolecule MoleculeGroup::viewAt(int idx) const
{
    idx = Index(idx).map(d->molviewidx_to_num.count());
    
    const tuple<MolNum,Index> &molviewidx = d->molviewidx_to_num.constData()[idx];
    
    PartialMolecule oldmol = d->molecules.at( molviewidx.get<0>(), molviewidx.get<1>() );

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the MolNum/Index index of the 'ith' view in this group

    \throw SireError::invalid_index
*/
const tuple<MolNum,Index>& MoleculeGroup::molViewIndexAt(int idx) const
{
    idx = Index(idx).map(d->molviewidx_to_num.count());
    
    return d->molviewidx_to_num.constData()[idx];
}

/** Return the number of the 'ith' molecule in this group

    \throw SireError::invalid_index
*/
MolNum MoleculeGroup::molNumAt(int idx) const
{
    return d->molidx_to_num.constData()[ Index(idx).map(d->molidx_to_num.count()) ];
}

/** Return the index of the view of the molecule viewed in 'molview'. This 
    is the index of this specific view, so you use this index with the
    MoleculeGroup::viewAt(int i) function. This returns -1 if this
    view is not in this group */
int MoleculeGroup::indexOf(const MoleculeView &molview) const
{
    //get this molecule
    const ViewsOfMol &mol = this->operator[](molview.data().number());
    
    //get the index of this view in the molecule
    int idx = mol.indexOf( molview.selection() );
    
    if (idx == -1)
        return idx;

    try
    {
        return d->molviewidx_to_num.indexOf( tuple<MolNum,Index>(mol.data().number(),
                                                                 Index(idx)) );
    }
    catch(...)
    {
        //the molecule layout has changed - it is definitely not in the group
        return -1;
    }
}

/** Return the index of the molecule with number 'molnum'. This is the index
    of the molecule itself, so you use this index with the MoleculeGroup::at(int i)
    function. This returns -1 if this molecule isn't in this group.
*/
int MoleculeGroup::indexOf(MolNum molnum) const
{
    return d->molidx_to_num.indexOf(molnum);
}

/** Return the views of the molecule with number 'molnum' from
    this group
    
    \throw SireMol::missing_molecule
*/
const ViewsOfMol& MoleculeGroup::molecule(MolNum molnum) const
{
    const ViewsOfMol &oldmol = d->molecules.at(molnum);
    
    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return the views of the molecule at index 'molidx'

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::molecule(MolIdx molidx) const
{
    return this->at( this->getMoleculeNumber(molidx) );
}

/** Return the views of the molecule with name 'molname'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
const ViewsOfMol& MoleculeGroup::molecule(const MolName &molname) const
{
    return this->at( this->getMoleculeNumber(molname) );
}

/** Return the views of the molecule that matches the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::molecule(const MolID &molid) const
{
    return this->at( this->getMoleculeNumber(molid) );
}

/** Return all of the molecules that match the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
Molecules MoleculeGroup::molecules(const MolID &molid) const
{
    QList<MolNum> molnums = this->map(molid);
    
    Molecules mols;
    
    foreach (MolNum molnum, molnums)
    {
        mols.add( this->at(molnum) );
    }
    
    return mols;
}

/** Obvious function used to shortcut the getMoleculeNumber(const MolID&)
    function
    
    \throw SireMol::missing_molecule
*/
MolNum MoleculeGroup::getMoleculeNumber(MolNum molnum) const
{
    this->assertContains(molnum);
    return molnum;
}

/** Return the number of the molecule at index 'molidx'

    \throw SireError::invalid_index
*/
MolNum MoleculeGroup::getMoleculeNumber(MolIdx molidx) const
{
    int i = molidx.map( d->molidx_to_num.count() );
    
    return d->molidx_to_num.constData()[i];
}

/** Return the number of the molecule with name 'molname'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
MolNum MoleculeGroup::getMoleculeNumber(const MolName &molname) const
{
    QList<MolNum> molnums = this->map(molname);
    
    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule with the name \"%1\"."
            "Their numbers are %2.")
                .arg(molname).arg(Sire::toString(molnums)),
                    CODELOC );

    return molnums.first();
}

/** Return the number of the molecule that matches the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
MolNum MoleculeGroup::getMoleculeNumber(const MolID &molid) const
{
    QList<MolNum> molnums = this->map(molid);
    
    if (molnums.count() > 1)
        throw SireMol::duplicate_molecule( QObject::tr(
            "There is more than one molecule that matches the ID %1. "
            "Matching molecules have numbers %2.")
                .arg(molid.toString())
                .arg(Sire::toString(molnums)), CODELOC );
                
    return molnums.first();
}

/** Obvious overload that shortcuts the map(const MolID&) function

    \throw SireMol::missing_molecule
*/
QList<MolNum> MoleculeGroup::map(MolNum molnum) const
{
    this->assertContains(molnum);
    
    QList<MolNum> molnums;
    molnums.append(molnum);
    
    return molnums;
}

/** Return the number of the molecule at index 'molidx'

    \throw SireError::invalid_index
*/
QList<MolNum> MoleculeGroup::map(MolIdx molidx) const
{
    int i = molidx.map( d->molidx_to_num.count() );

    QList<MolNum> molnums;
    molnums.append( d->molidx_to_num.constData()[i] );
    
    return molnums;
}

/** Return the numbers of the molecules that are called 'molname'

    \throw SireMol::missing_molecule
*/
QList<MolNum> MoleculeGroup::map(const MolName &molname) const
{
    return molname.map( this->molecules() );
}

/** Return the numbers of the molecules that match the ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> MoleculeGroup::map(const MolID &molid) const
{
    return molid.map(*this);
}

/** Return whether or not this group contains any views of the 
    molecule with number 'molnum' */
bool MoleculeGroup::contains(MolNum molnum) const
{
    return d->molecules.contains(molnum);
}

/** Return whether or not this group contains a molecule at index 'molidx' */
bool MoleculeGroup::contains(MolIdx molidx) const
{
    try
    {
        molidx.map( d->molidx_to_num.count() );
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not this group contains any molecules called 'molname' */
bool MoleculeGroup::contains(const MolName &molname) const
{
    for (Molecules::const_iterator it = d->molecules.constBegin();
         it != d->molecules.constEnd();
         ++it)
    {
        if (it->data().name() == molname)
            return true;
    }

    return false;
}

/** Return whether or not this group contains any molecules that
    match the ID 'molid' */
bool MoleculeGroup::contains(const MolID &molid) const
{
    try
    {
        molid.map(*this);
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not this group contains any version of 
    the view of the molecule in 'molview' */
bool MoleculeGroup::contains(const MoleculeView &molview) const
{
    return d->molecules.contains(molview);
}

/** Return whether or not this group contains all of the views
    of any version of the molecule in 'molviews' */
bool MoleculeGroup::contains(const ViewsOfMol &molviews) const
{
    return d->molecules.contains(molviews);
}

/** Return whether or not this group contains all of the 
    views of any version of all of the molecules contained
    in 'molecules' */
bool MoleculeGroup::contains(const Molecules &molecules) const
{
    return d->molecules.contains(molecules);
}

/** Return whether or not this group contains all of the 
    views of any version of all of the molecules contained
    in the group 'other' */
bool MoleculeGroup::contains(const MoleculeGroup &other) const
{
    return this->contains(other.molecules());
}

/** Return whether or not this group contains any version 
    of any of the atoms of the molecule in 'molview' */
bool MoleculeGroup::intersects(const MoleculeView &molview) const
{
    return d->molecules.intersects(molview);
}

/** Return whether or not this group contains any version
    of any of the atoms in any of the molecules in 'molecules' */
bool MoleculeGroup::intersects(const Molecules &other) const
{
    return d->molecules.intersects(other);
}

/** Return whether or not this group contains any version
    of any of the atoms in any of the molecules contained in 
    the group 'other' */
bool MoleculeGroup::intersects(const MoleculeGroup &other) const
{
    return this->intersects(other.molecules());
}

/** Return the number of molecules in this group */
int MoleculeGroup::nMolecules() const
{
    return d->molidx_to_num.count();
}

/** Return the number of views of molecules in this group - 
    this must always be greater or equal to the number of 
    molecules! */
int MoleculeGroup::nViews() const
{
    return d->molviewidx_to_num.count();
}

/** Return the number of views of the molecule with number 'molnum'
    that are present in this group.
    
    \throw SireMol::missing_molecule
*/
int MoleculeGroup::nViews(MolNum molnum) const
{
    return d->molecules.nViews(molnum);
}

/** Return the number of views of the molecule(s) that match
    the ID 'molid'
    
    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
int MoleculeGroup::nViews(const MolID &molid) const
{
    QList<MolNum> molnums = this->map(molid);
    
    int nviews = 0;
    
    foreach (MolNum molnum, molnums)
    {
        nviews += this->nViews(molnum);
    }
    
    return nviews;
}

/** Return the number of views of the molecule at index 'idx'
    in this group */
int MoleculeGroup::nViews(Index idx) const
{
    return this->nViews( d->molidx_to_num.constData()
                             [idx.map(d->molidx_to_num.count())] );
}

/** Return whether or not this group is empty */
bool MoleculeGroup::isEmpty() const
{
    return d->molecules.isEmpty();
}

/** Return all views of all of the molecules in this group */
const Molecules& MoleculeGroup::molecules() const
{
    if (this->needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules;
}

/** Return a reference to the first molecule in the group 

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::first() const
{
    const ViewsOfMol &oldmol = d->molecules.first();

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return a reference to the last molecule in the group

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::last() const
{
    const ViewsOfMol &oldmol = d->molecules.last();

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return a reference to the first molecule in the group

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::front() const
{
    const ViewsOfMol &oldmol = d->molecules.front();

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return a reference to the last molecule in the group 

    \throw SireError::invalid_index
*/
const ViewsOfMol& MoleculeGroup::back() const
{
    const ViewsOfMol &oldmol = d->molecules.back();

    if (workspace.isEmpty())
        return oldmol;
    else
        return workspace.getUpdated(oldmol);
}

/** Return an iterator pointing to the first molecule
    in the group */
MoleculeGroup::const_iterator MoleculeGroup::begin() const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules.begin();
}

/** Return an iterator pointing one space after the last
    molecule in the group */
MoleculeGroup::const_iterator MoleculeGroup::end() const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules.end();
}

/** Return an iterator pointing to the first molecule
    in the group */
MoleculeGroup::const_iterator MoleculeGroup::constBegin() const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules.constBegin();
}

/** Return an iterator pointing to one space after
    the last molecule in the group */
MoleculeGroup::const_iterator MoleculeGroup::constEnd() const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules.constEnd();
}

/** Return an iterator pointing to the molecule with 
    number 'molnum'. If there is no such molecule then
    this returns MoleculeGroup::end() */
MoleculeGroup::const_iterator MoleculeGroup::find(MolNum molnum) const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules.find(molnum);
}

/** Return an iterator pointing to the molecule with 
    number 'molnum'. If there is no such molecule then
    this returns MoleculeGroup::end() */
MoleculeGroup::const_iterator MoleculeGroup::constFind(MolNum molnum) const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return d->molecules.constFind(molnum);
}

/** Return an iterator that points to the molecule that matches
    the ID 'molid'
    
    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
MoleculeGroup::const_iterator MoleculeGroup::find(const MolID &molid) const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return this->find( this->getMoleculeNumber(molid) );
}

/** Return an iterator that points to the molecule that matches
    the ID 'molid'
    
    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
MoleculeGroup::const_iterator MoleculeGroup::constFind(const MolID &molid) const
{
    if (needsAccepting())
    {
        const_cast<MoleculeGroup*>(this)->accept();
    }

    return this->constFind( this->getMoleculeNumber(molid) );
}

/** Return the numbers of all molecules present in this group,
    in the order that the molecules appear in this group */
const QVector<MolNum>& MoleculeGroup::molNums() const
{
    return d->molidx_to_num;
}

/** Return the numbers and view indicies of all views in this group,
    in the order that the views appear in this group */
const QVector< tuple<MolNum,Index> >& MoleculeGroup::molViewIndicies() const
{
    return d->molviewidx_to_num;
}

/** Return the set of all names of the molecules in this group */
QSet<MolName> MoleculeGroup::molNames() const
{
    QSet<MolName> molnames;
    
    for (Molecules::const_iterator it = d->molecules.constBegin();
         it != d->molecules.constEnd();
         ++it)
    {
        molnames.insert( it->data().name() );
    }

    return molnames;
}

/** Assert that this group contains a view of any part of the 
    molecule with number 'molnum'
    
    \throw SireMol::missing_molecule
*/
void MoleculeGroup::assertContains(MolNum molnum) const
{
    d->molecules.assertContains(molnum);
}

/** Assert that this group contains a molecule called 'molname'

    \throw SireMol::missing_molecule
*/
void MoleculeGroup::assertContains(const MolName &molname) const
{
    if (not this->contains(molname))
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule called \"%1\" in this group. Use "
            "MoleculeGroup::molNames() to get the set of names of molecules "
            "that are in this group.")
                .arg(molname), CODELOC );
}

/** Return the name of this group */
const MGName& MoleculeGroup::name() const
{
    return d->name;
}

/** Return the ID number of this group */
MGNum MoleculeGroup::number() const
{
    return d->number;
}

/** Return the version number of this group */
const Version& MoleculeGroup::version() const
{
    if (workspace.isEmpty())
        return d->version.version();
    else
        return workspace.version().version();
}

/** Change the name of this group */
void MoleculeGroup::setName(const QString &new_name)
{
    if (new_name != this->name())
    {
        d->name = MGName(new_name);
        d->incrementMajor();
    }
}

/** Change the number of this group */
void MoleculeGroup::setNumber(quint32 new_number)
{
    MGNum new_num(new_number);

    if (new_num != this->number())
    {
        d->number = new_num;
        d->version = mgnum_registry.registerObject(d->number);
    }
}

/** Give this group a new, unique number */
void MoleculeGroup::setNewNumber()
{
    MGNum new_num = MGNum::getUniqueNumber();
    this->setNumber(new_num);
}

/** Return the major version number of this group. This number
    changes whenever views are added or removed from this group,
    or when the name of this group changes */
quint64 MoleculeGroup::majorVersion() const
{
    return this->version().majorVersion();
}

/** Return the minor version number of this group. This number
    changes whenever any of the versions of molecules in this group
    are changed. This number is reset to zero whenever the major
    version number of this group is changed. */
quint64 MoleculeGroup::minorVersion() const
{
    return this->version().minorVersion();
}

/** Add the view of the molecule in 'molview' to this group. 
    This adds the view as a duplicate if it already exists 
    in this group */
void MoleculeGroup::add(const MoleculeView &molview)
{
    if (molview.selection().isEmpty())
        return;

    accept();

    MolNum molnum = molview.data().number();

    MolGroupPvt &dref = *d;

    if (not dref.molecules.contains(molnum))
        dref.molidx_to_num.append(molnum);
    
    dref.molecules.add(molview);
   
    dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,
                                         Index(dref.molecules.nViews(molnum) - 1)) );

    dref.incrementMajor();
}

/** Add the views of the molecule in 'molviews' to this group.
    This adds the views as duplicates if they already exist
    in this group */
void MoleculeGroup::add(const ViewsOfMol &molviews)
{
    if (molviews.isEmpty())
        return;
    else if (molviews.nViews() == 1)
    {
        this->add(molviews.at(0));
        return;
    }

    accept();

    MolNum molnum = molviews.number();
    
    MolGroupPvt &dref = *d;
    
    if (not dref.molecules.contains(molnum))
        dref.molidx_to_num.append(molnum);
        
    dref.molecules.add(molviews);
    
    quint32 nviews = dref.molecules.nViews(molnum);
    
    for (quint32 i = nviews - molviews.count(); i < nviews; ++i)
    {
        dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,Index(i)) );
    }
    
    dref.incrementMajor();
}

/** Add all of the molecules in 'molecules' to this group.
    This duplicates any molecules that already exist in this 
    group. */
void MoleculeGroup::add(const Molecules &molecules)
{
    if (molecules.isEmpty())
        return;
    else if (molecules.count() == 1)
    {
        this->add(molecules.first());
        return;
    }

    accept();

    MolGroupPvt &dref = *d;
    
    //add the molecules to the index
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        MolNum molnum = it.key();
        
        quint32 nviews = 0;
        
        if (not dref.molecules.contains(molnum))
            dref.molidx_to_num.append(molnum);
        else
            nviews = dref.molecules.nViews(molnum);

        quint32 n_newviews = it->nViews();
        
        for (quint32 i=nviews; i<nviews + n_newviews; ++i)
        {
            dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,Index(i)) );
        }
    }

    //now add the molecules themselves
    dref.molecules.add(molecules);
    
    dref.incrementMajor();
}

/** Add the molecules in 'MoleculeGroup' to this set. This adds the 
    molecules and views in the same order as they appear in 
    'MoleculeGroup', adding them as duplicates if they already
    exist in this set. Note that the version of the molecule
    will be taken from this set. */
void MoleculeGroup::add(const MoleculeGroup &molgroup)
{
    if (molgroup.isEmpty())
        return;
    
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        copy.accept();
        this->add(copy);
        return;
    }
    
    if (this->isEmpty())
    {
        MolGroupPvt &dref = *d;
        dref.molecules = molgroup.d->molecules;
        dref.molidx_to_num = molgroup.d->molidx_to_num;
        dref.molviewidx_to_num = molgroup.d->molviewidx_to_num;
        
        dref.incrementMajor();
        return;
    }
    
    accept();
    
    SharedDataPointer<MolGroupPvt> old_state = d;
    
    try
    {
        MolGroupPvt &dref = *d;
    
        //append the other group's index onto this group
        foreach (MolNum molnum, molgroup.d->molidx_to_num)
        {
            if (not dref.molecules.contains(molnum))
                dref.molidx_to_num.append(molnum);
        }
    
        for (QVector< tuple<MolNum,Index> >::const_iterator 
                                it = molgroup.d->molviewidx_to_num.constBegin();
             it != molgroup.d->molviewidx_to_num.constEnd();
             ++it)
        {
            MolNum molnum = it->get<0>();

            int nviews = 0;
            
            if (dref.molecules.contains(molnum))
                nviews = dref.molecules.nViews(molnum);
    
            dref.molviewidx_to_num.append( tuple<MolNum,Index>( molnum,
                            Index(it->get<1>() + nviews) ) );
        }
    
        //now add the molecules themselves to this set
        dref.molecules.add(molgroup.d->molecules);
    
        dref.incrementMajor();
    }
    catch(...)
    {
        d = old_state;
        throw;
    }
}

/** Add the view of the molecule in 'molview' to this group. 
    This only adds the view if it does not already exist in 
    this group, and returns whether or not the view was added */
bool MoleculeGroup::addIfUnique(const MoleculeView &molview)
{
    if (molview.selection().isEmpty())
        return false;

    accept();

    MolNum molnum = molview.data().number();

    MolGroupPvt &dref = *d;
    
    bool hasmol = dref.molecules.contains(molnum);
    
    if (dref.molecules.addIfUnique(molview))
    {
        //the view has been added!
        if (not hasmol)
            dref.molidx_to_num.append(molnum);
            
        dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,
                                         Index(dref.molecules.nViews(molnum) - 1)) );
    
        dref.incrementMajor();
        return true;
    }
    else
        return false;
}

/** Add the views of the molecule 'molviews' to this group,
    adding the views only if they don't already exist in this
    group. This returns the views that were added successfully
    to this group. */
ViewsOfMol MoleculeGroup::addIfUnique(const ViewsOfMol &molviews)
{
    if (molviews.isEmpty())
        return ViewsOfMol();
    
    accept();
    
    MolNum molnum = molviews.number();
    
    MolGroupPvt &dref = *d;
    
    bool hasmol = dref.molecules.contains(molnum);
    
    ViewsOfMol added_views = dref.molecules.addIfUnique(molviews);
    
    if (not added_views.isEmpty())
    {
        if (not hasmol)
            dref.molidx_to_num.append(molnum);
            
        quint32 nviews = dref.molecules.nViews(molnum);
        quint32 nadded = added_views.nViews();
        
        for (quint32 i=nviews-nadded-1; i<nviews; ++i)
        {
            dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,Index(i)) );
        }
        
        dref.incrementMajor();
    }
    
    return added_views;
}

/** Add the views of the molecules in 'molecules' to this group. This
    only adds views that don't already exist in this group. This
    returns all of the views that were successfully added. */
QList<ViewsOfMol> MoleculeGroup::addIfUnique(const Molecules &molecules)
{
    if (molecules.isEmpty())
        return QList<ViewsOfMol>();
    
    accept();
    
    MolGroupPvt &dref = *d;

    //which molecules already exist in this group?
    QSet<MolNum> groupmols = dref.molecules.molNums();
    
    //add the molecules to this group
    QList<ViewsOfMol> added_mols = dref.molecules.addIfUnique(molecules);
    
    //now update the index
    if (not added_mols.isEmpty())
    {
        foreach (const ViewsOfMol &added_mol, added_mols)
        {
            MolNum molnum = added_mol.number();
            
            if (not groupmols.contains(molnum))
                dref.molidx_to_num.append(molnum);
                
            quint32 nviews = dref.molecules.nViews(molnum);
            quint32 nadded = added_mol.nViews();
            
            for (quint32 i=nviews-nadded-1; i<nviews; ++i)
            {
                dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,Index(i)) );
            }
        }
        
        dref.incrementMajor();
    }

    return added_mols;
}

/** Add the views/molecules in 'MoleculeGroup' to this group, but
    only if they don't already exist in this group. This has
    the same action as MoleculeGroup::addIfUnique(molecules), but
    it ensures that the added views are in the same order as
    in 'MoleculeGroup'. This is costly, so if you don't care
    about the added order, then use 
    MoleculeGroup::addIfUnique(MoleculeGroup.molecules()) instead. 
    
    This returns the added views.
*/
QList<ViewsOfMol> MoleculeGroup::addIfUnique(const MoleculeGroup &molgroup)
{
    QList<ViewsOfMol> added_mols;

    if (molgroup.isEmpty())
        return added_mols;
    
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        copy.accept();
        return this->addIfUnique(copy);
    }
    
    accept();
    
    MolGroupPvt &dref = *d;

    //which molecules already exist in this group?
    QSet<MolNum> groupmols = dref.molecules.molNums();
    
    for (QVector< tuple<MolNum,Index> >::const_iterator 
                            it = dref.molviewidx_to_num.constBegin();
         it != dref.molviewidx_to_num.constEnd();
         ++it)
    {
        MolNum molnum = it->get<0>();
        Index i = it->get<1>();
    
        PartialMolecule molview = molgroup[molnum][i];
        
        if (dref.molecules.addIfUnique(molview))
            added_mols.append(molview);
    }

    //now update the index
    if (not added_mols.isEmpty())
    {
        foreach (const ViewsOfMol &added_mol, added_mols)
        {
            MolNum molnum = added_mol.number();
            
            if (not groupmols.contains(molnum))
                dref.molidx_to_num.append(molnum);
                
            quint32 nviews = dref.molecules.nViews(molnum);
            quint32 nadded = added_mol.nViews();
            
            for (quint32 i=nviews-nadded-1; i<nviews; ++i)
            {
                dref.molviewidx_to_num.append( tuple<MolNum,Index>(molnum,Index(i)) );
            }
        }
        
        dref.incrementMajor();
    }
        
    return added_mols;
}

/** Synonym for MoleculeGroup::addIfUnique(molview) */
bool MoleculeGroup::unite(const MoleculeView &molview)
{
    return this->addIfUnique(molview);
}

/** Synonym for MoleculeGroup::addIfUnique(molviews) */
ViewsOfMol MoleculeGroup::unite(const ViewsOfMol &molviews)
{
    return this->addIfUnique(molviews);
}

/** Synonym for MoleculeGroup::addIfUnique(molecules) */
QList<ViewsOfMol> MoleculeGroup::unite(const Molecules &molecules)
{
    return this->addIfUnique(molecules);
}

/** Synonym for MoleculeGroup::addIfUnique(MoleculeGroup). The 
    function MoleculeGroup::addIfUnique(MoleculeGroup.molecules()) is
    quicker if you don't care about the order in which
    the views are added. */
QList<ViewsOfMol> MoleculeGroup::unite(const MoleculeGroup &MoleculeGroup)
{
    return this->addIfUnique(MoleculeGroup);
}

bool MoleculeGroup::_pvt_remove(const MoleculeView &molview)
{
    AtomSelection selected_atoms = molview.selection();

    if (selected_atoms.isEmpty())
        return false;

    accept();

    MolNum molnum = molview.data().number();

    Molecules::const_iterator it = d.constData()->molecules.find(molnum);

    if (it == d.constData()->molecules.end())
        return false;
        
    int viewidx = it->indexOf(selected_atoms);
    
    if (viewidx == -1)
        //this view is not present in this group
        return false;
    
    //the view is present - remove it!    
    MolGroupPvt &dref = *d;
    
    //remove all of the views of this molecule...
    ViewsOfMol molviews = dref.molecules.remove(molview.data().number());
    
    //remove the specified view from this set...
    molviews.removeAt(viewidx);
    
    //if there are any views left then add them back
    //to the collection of molecules
    if (not molviews.isEmpty())
    {
        dref.molecules.add(molviews);
        
        //the molecule has only been partially removed,
        //so remove this view from the molview index
        QMutableVectorIterator< tuple<MolNum,Index> > it(dref.molviewidx_to_num);
        
        while (it.hasNext())
        {
            tuple<MolNum,Index> &molviewidx = it.next();
            
            if (molviewidx.get<0>() == molnum)
            {
                if (molviewidx.get<1>() == viewidx)
                    it.remove();
                else if (molviewidx.get<1>() > viewidx)
                    molviewidx = tuple<MolNum,Index>(molnum,
                                                     Index(molviewidx.get<1>()-1));
            }
        }
    }
    else
    {
        //this molecule has been completely removed...
        //remove it completely from the index
        dref.molidx_to_num.remove( dref.molidx_to_num.indexOf(molnum) );
        
        QMutableVectorIterator< tuple<MolNum,Index> > it(dref.molviewidx_to_num);
        
        while (it.hasNext())
        {
            const tuple<MolNum,Index> &molviewidx = it.next();
        
            if (molviewidx.get<0>() == molnum)
                it.remove();
        }
    }
    
    return true;
}

/** Remove the view of the molecule in 'molview' from this set.
    This only removes the first such view from the set, and 
    returns whether or not any view was removed */
bool MoleculeGroup::remove(const MoleculeView &molview)
{
    if (this->_pvt_remove(molview))
    {
        d->incrementMajor();
        return true;
    }
    else
        return false;
}

ViewsOfMol MoleculeGroup::_pvt_remove(const ViewsOfMol &molviews)
{
    if (molviews.isEmpty())
        return molviews;
    
    accept();
    
    int nviews = molviews.nViews();
    
    QList<AtomSelection> removed_views;
    
    for (int i=0; i<nviews; ++i)
    {
        PartialMolecule view = molviews.at(i);
    
        if (this->_pvt_remove(view))
        {
            removed_views.append(view.selection());
        }
    }
    
    if (removed_views.isEmpty())
        return ViewsOfMol();
    else
        return ViewsOfMol(molviews.data(), removed_views);
}

/** Remove all of the views of the molecule in 'molviews' from this
    set. This only removes the first such view of any duplicates
    from this set, and returns the views that were removed */
ViewsOfMol MoleculeGroup::remove(const ViewsOfMol &molviews)
{
    ViewsOfMol removed_views = this->_pvt_remove(molviews);
    
    if (not removed_views.isEmpty())
        d->incrementMajor();
        
    return removed_views;
}

/** Remove all of the molecules listed in 'molecules' from this set. 
    This only removes the first of any duplicated views in this set.
    This returns the views/molecules that were successfully removed. */
QList<ViewsOfMol> MoleculeGroup::remove(const Molecules &molecules)
{
    QList<ViewsOfMol> removed_mols;
    
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        ViewsOfMol removed_views = this->_pvt_remove(*it);
        
        if (not removed_views.isEmpty())
            removed_mols.append(removed_views);
    }

    if (not removed_mols.isEmpty())
        d->incrementMajor();
    
    return removed_mols;
}

/** Remove all of the molecules from the group 'MoleculeGroup' from this set.
    This only removes the first of any duplicated views in this set.
    This returns the views/molecules that were sucessfully removed. */
QList<ViewsOfMol> MoleculeGroup::remove(const MoleculeGroup &molgroup)
{
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        copy.accept();
        return this->remove(copy);
    }

    return this->remove(molgroup.molecules());
}

bool MoleculeGroup::_pvt_removeAll(const MoleculeView &molview)
{
    accept();

    //just keep removing it until the view has gone completely!
    bool removed_a_view = false;
    
    while (this->_pvt_remove(molview))
    {
        removed_a_view = true;
    }
    
    return removed_a_view;
}

/** Remove all copies of the view of the molecule in 'molview' from this
    group. This removes all copies if this view is duplicated in this
    group, and returns whether or not any views were removed. */
bool MoleculeGroup::removeAll(const MoleculeView &molview)
{
    if (this->_pvt_removeAll(molview))
    {
        d->incrementMajor();
        return true;
    }
    else
        return false;
}

ViewsOfMol MoleculeGroup::_pvt_removeAll(const ViewsOfMol &molviews)
{
    if (molviews.isEmpty())
        return molviews;
    
    accept();
    
    QList<AtomSelection> removed_views;
    
    int nviews = molviews.nViews();
    
    for (int i=0; i<nviews; ++i)
    {
        PartialMolecule view = molviews.at(i);
        
        if (this->_pvt_removeAll(view))
            removed_views.append(view.selection());
    }
    
    if (removed_views.isEmpty())
        return ViewsOfMol();
    else
        return ViewsOfMol(molviews.data(), removed_views);
}

/** Remove all copies of all of the views of the molecule in 'molviews'.
    This removes all copies of any duplicated views in this group,
    and returns the views that were successfully removed. */
ViewsOfMol MoleculeGroup::removeAll(const ViewsOfMol &molviews)
{
    ViewsOfMol removed_views = this->_pvt_remove(molviews);
    
    if (not removed_views.isEmpty())
        d->incrementMajor();
        
    return removed_views;
}

/** Remove all copies of all of the views of the molecules in 'molecules'.
    This removes all copies of any duplicated views in this group.
    This returns the molecules/views that were removed. */
QList<ViewsOfMol> MoleculeGroup::removeAll(const Molecules &molecules)
{
    QList<ViewsOfMol> removed_mols;
    
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        ViewsOfMol removed_views = this->_pvt_removeAll(*it);
        
        if (not removed_views.isEmpty())
            removed_mols.append(removed_views);
    }
    
    if (not removed_mols.isEmpty())
        d->incrementMajor();
    
    return removed_mols;
}

/** Remove all copies of all of the views of the molecules in the
    group 'MoleculeGroup'. This removes all copies of any duplicated 
    views in this group. This returns the molecules/views that
    were removed
*/
QList<ViewsOfMol> MoleculeGroup::removeAll(const MoleculeGroup &molgroup)
{
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        copy.accept();
        return removeAll(copy);
    }

    return this->removeAll(molgroup.molecules());
}

ViewsOfMol MoleculeGroup::_pvt_remove(MolNum molnum)
{
    if (not this->contains(molnum))
        return ViewsOfMol();

    accept();

    MolGroupPvt &dref = *d;

    //remove the molecule
    ViewsOfMol removed_views = dref.molecules.remove(molnum);
    
    //now remove it from the index
    dref.molidx_to_num.remove( dref.molidx_to_num.indexOf(molnum) );
    
    QMutableVectorIterator< tuple<MolNum,Index> > it(dref.molviewidx_to_num);
    
    while (it.hasNext())
    {
        const tuple<MolNum,Index> &molviewidx = it.next();
        
        if (molviewidx.get<0>() == molnum)
            it.remove();
    }
    
    return removed_views;
}

/** Completely remove all views of the molecule with number 'molnum'
    from this group. This returns the views that were removed */
ViewsOfMol MoleculeGroup::remove(MolNum molnum)
{
    ViewsOfMol removed_views = this->_pvt_remove(molnum);
    
    if (not removed_views.isEmpty())
        d->incrementMajor();
        
    return removed_views;
}

/** Remove all views of the molecules whose numbers are in 'molnums'.
    This returns the views that were removed. */
QList<ViewsOfMol> MoleculeGroup::remove(const QSet<MolNum> &molnums)
{
    QList<ViewsOfMol> removed_mols;

    foreach (MolNum molnum, molnums)
    {
        ViewsOfMol removed_views = this->_pvt_remove(molnum);
        
        if (not removed_views.isEmpty())
            removed_mols.append(removed_views);
    }
    
    if (not removed_mols.isEmpty())
        d->incrementMajor();
    
    return removed_mols;
}

/** Remove all of the molecules from this group */
void MoleculeGroup::removeAll()
{
    if (not this->isEmpty())
    {
        MolGroupPvt &dref = *d;
        
        dref.molecules.clear();
        dref.molidx_to_num.clear();
        dref.molviewidx_to_num.clear();
        
        workspace.clear();
        
        dref.incrementMajor();
    }
}

/** Update this group so that the molecule in this group whose 
    data is in 'moldata' is also at the same version as 'moldata'.
    
    This does nothing if there is no such molecule in this 
    group, or if it is already at this version, and this returns
    whether or not this changes the group. */
bool MoleculeGroup::update(const MoleculeData &moldata, bool auto_commit)
{
    if (auto_commit)
    {
        if (this->needsAccepting())
            this->accept();
        
        if (d.constData()->molecules.contains(moldata.number()))
        {
            if (d.constData()->molecules.at(moldata.number()).version() != moldata.version())
            {
                d->molecules.update(moldata);
                d->version.incrementMinor();
                return true;
            }
        }

        return false;
    }
    else
    {
        if (d.constData()->molecules.contains(moldata.number()))
        {
            if (d.constData()->molecules.at(moldata.number()).version() != moldata.version())
            {
                if (workspace.isEmpty())
                {
                    workspace.setVersion(d.constData()->version);
                }
                
                workspace.push(moldata);
                workspace.incrementMinor();
                return true;
            }
        }
        return false;
    }
}

/** Update this group so that the molecule in this group that
    is also viewed in 'molview' is updated to the same
    version as 'molview'.
    
    This does nothing if there is no such molecule in this 
    group, or if it is already at this version, and this returns
    whether or not this changes the group. */
bool MoleculeGroup::update(const MoleculeView &molview, bool auto_commit)
{
    return this->update(molview.data(), auto_commit);
}

/** Update this group so that the contained molecules have the 
    same versions as the molecules in 'molecules'. This does
    nothing if none of these molecules are in this group, or
    if they are already at the same versions. This returns
    the list of molecules that were changed by this update. */
QList<Molecule> MoleculeGroup::update(const Molecules &molecules, bool auto_commit)
{
    QList<Molecule> updated_mols;

    if (auto_commit)
    {
        if (this->needsAccepting())
            this->accept();
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (d.constData()->molecules.contains(it.key()))
            {
                if (d.constData()->molecules.at(it.key()).version() != it.value().version())
                {
                    d->molecules.update(it.value().data());
                    updated_mols.append(it.value().molecule());
                }
            }
        }
        
        if (not updated_mols.isEmpty())
            d->version.incrementMinor();
    }
    else
    {
        bool must_create_version = false;
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (d.constData()->molecules.contains(it.key()))
            {
                if (d.constData()->molecules.at(it.key()).version() != it.value().version())
                {
                    if (workspace.isEmpty())
                        must_create_version = true;
                
                    workspace.push(it.value().data());
                    updated_mols.append(it.value().molecule());
                }
            }
        }
        
        if (not updated_mols.isEmpty())
        {
            if (must_create_version)
                workspace.setVersion( d.constData()->version );
            
            workspace.incrementMinor();
        }
    }
    
    return updated_mols;
}

/** Update this group so that it has the same version of molecules
    as those in 'MoleculeGroup'. This does nothing if this group
    doesn't contain any of the molecules in 'MoleculeGroup', or
    if it already has the molecules at the same version. This
    returns the list of molecules that were changed by this
    update
*/
QList<Molecule> MoleculeGroup::update(const MoleculeGroup &molgroup, bool auto_commit)
{
    if (molgroup == *this)
        //there is nothing to update!
        return QList<Molecule>();
    else
    {
        if (molgroup.needsAccepting())
        {
            MoleculeGroup copy(molgroup);
            copy.accept();
            return this->update(copy.molecules(), auto_commit);
        }
        else
            return this->update(molgroup.molecules(), auto_commit);
    }
}

void MoleculeGroup::_pvt_setContents(const Molecules &molecules)
{
    workspace.clear();

    MolGroupPvt &dref = *d;
    
    //set the molecules
    dref.molecules = molecules;
    
    //resize the two indexes
    dref.molidx_to_num.resize(molecules.nMolecules());
    dref.molviewidx_to_num.resize(molecules.nViews());

    MolNum *molidx_to_num_array = dref.molidx_to_num.data();
    tuple<MolNum,Index> *molviewidx_to_num_array = dref.molviewidx_to_num.data();
    
    int imol = 0;
    int iview = 0;
    
    //index the molecules and views...
    for (Molecules::const_iterator it = molecules.begin();
         it != molecules.end();
         ++it)
    {
        const ViewsOfMol &molviews = *it;
        
        //index the molecule
        molidx_to_num_array[imol] = it.key();
        ++imol;
        
        //now index its views
        int nviews = molviews.nViews();
        
        for (int i=0; i < nviews; ++i)
        {
            molviewidx_to_num_array[iview+i] = tuple<MolNum,Index>(it.key(), Index(i));
        }
        
        iview += nviews;
    }
    
    dref.incrementMajor();
}

/** Set the contents of this group to 'molecules'. This clears
    any existing contents of this group. */
bool MoleculeGroup::setContents(const Molecules &molecules)
{
    workspace.clear();

    if (molecules.isEmpty())
    {
        if (this->isEmpty())
            return false;
        else
        {
            this->removeAll();
            return true;
        }
    }
    else if (this->isEmpty())
    {
        this->_pvt_setContents(molecules);
        return true;
    }

    const Molecules &oldmols = this->molecules();

    if (oldmols == molecules)
        //there is nothing to do
        return false;
        
    //see if the actual list of molecules has changed, or whether
    //it is just their versions...
    if (oldmols.nMolecules() == molecules.nMolecules())
    {
        bool different = false;
    
        //loop over all molecules and views and check that the same
        //molecules are in both sets, with the same views
        for (Molecules::const_iterator it = oldmols.begin();
             it != oldmols.end();
             ++it)
        {
            if (not molecules.contains(it.key()))
            {
                different = true;
                break;
            }
        
            const ViewsOfMol &oldviews = *it;
            const ViewsOfMol &newviews = molecules.molecule(it.key());
            
            int nviews = oldviews.nViews();
            
            if (nviews != newviews.nViews())
            {
                different = true;
                break;
            }
            
            for (int i=0; i<nviews; ++i)
            {
                if (oldviews.selection(i) != newviews.selection(i))
                {
                    different = true;
                    break;
                }
            }
            
            if (different)
                break;
        }
        
        if (not different)
        {
            //ok - the same molecules and views are in the two sets.
            //All that's changed are the versions :-)
            //Just update the molecules and increment the minor version 
            //counter
            d->molecules = molecules;
            d->incrementMinor();
            return true;
        }
    }
    
    //ok, the two sets of molecules contain different molecules and/or views.
    //We'll just have to reindex them completely.
    this->_pvt_setContents(molecules);
    return true;
}

/** Set the contents of this group so that it only contains the 
    view 'molview'. This clears any existing contents of this group */
bool MoleculeGroup::setContents(const MoleculeView &molview)
{
    return this->setContents( Molecules(molview) );
}

/** Set the contents of this group so that it only contains the 
    views of the molecule in 'molviews'. This clears any existing
    contents of this group. */
bool MoleculeGroup::setContents(const ViewsOfMol &molviews)
{
    return this->setContents( Molecules(molviews) );
}

/** Set the contents of this group so that it equals that
    of the group 'MoleculeGroup'. This sets the contents and
    also preserves the same order of molecules/views as
    in 'MoleculeGroup'
*/
bool MoleculeGroup::setContents(const MoleculeGroup &molgroup)
{
    if (molgroup.needsAccepting())
    {
        MoleculeGroup copy(molgroup);
        return this->setContents(copy);
    }

    //if the other group has the same ID number then we are reverting
    //this group back to another version
    if (this->number() == molgroup.number())
    {
        bool changed = (this->majorVersion() != molgroup.majorVersion() or
                        this->minorVersion() != molgroup.minorVersion());

        d = molgroup.d;
        
        return changed;
    }
        
    if (d.constData()->molecules == molgroup.d->molecules and
        d.constData()->molidx_to_num == molgroup.d->molidx_to_num and
        d.constData()->molviewidx_to_num == molgroup.d->molviewidx_to_num)
    {
        //the contents are exactly the same!
        return false;
    }
    
    MolGroupPvt &dref = *d;
    
    dref.molecules = molgroup.d->molecules;
    dref.molidx_to_num = molgroup.d->molidx_to_num;
    dref.molviewidx_to_num = molgroup.d->molviewidx_to_num;
    
    dref.incrementMajor();
    
    return true;
}

/** Return whether or not this molecule group has a temporary workspace that needs accepting */
bool MoleculeGroup::needsAccepting() const
{
    return not workspace.isEmpty();
}

/** Tell the molecule group that the last move was accepted. This tells the
    group to make permanent any temporary changes that were used a workspace
    to avoid memory allocation during a move */
void MoleculeGroup::accept()
{
    if (needsAccepting())
    {
        for (int i=0; i<workspace.count(); ++i)
        {
            d->molecules.update( workspace.constData()[i].read() );
        }
        
        d->version = workspace.version();
        workspace.clear();
    }
}

const char* MoleculeGroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MoleculeGroup>() );
}

MoleculeGroup* MoleculeGroup::clone() const
{
    return new MoleculeGroup(*this);
}
