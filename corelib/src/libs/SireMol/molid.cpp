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

#include "molid.h"
#include "molidx.h"
#include "molnum.h"
#include "molname.h"

#include "molatomid.h"

#include "atomidx.h"

#include "specifymol.h"

#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"
#include "selectormol.h"

#include "mgidx.h"

#include "mover.hpp"

#include "SireBase/incremint.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

///////
/////// Implementation of MolID
///////

/** Constructor */
MolID::MolID() : SireID::ID()
{}

/** Copy constructor */
MolID::MolID(const MolID &other) : SireID::ID(other)
{}

/** Destructor */
MolID::~MolID()
{}

/** Specify the ith molecule that matches this ID */
SpecifyMol MolID::operator[](int i) const
{
    return SpecifyMol(*this, i);
}

/** Specify the ith molecule that matches this ID */
SpecifyMol MolID::operator()(int i) const
{
    return SpecifyMol(*this, i);
}

/** Specify the ith to jth molecules that match this ID */
SpecifyMol MolID::operator()(int i, int j) const
{
    return SpecifyMol(*this, i, j);
}

/** Combine this ID with another molecule ID */
IDAndSet<MolID> MolID::operator+(const MolID &other) const
{
    return IDAndSet<MolID>(*this, other);
}

/** Syntactic sugar for operator+ */
IDAndSet<MolID> MolID::operator&&(const MolID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
IDAndSet<MolID> MolID::operator&(const MolID &other) const
{
    return this->operator+(other);
}

/** Search for matching atoms in the matching molecules */
MolAtomID MolID::operator+(const AtomID &other) const
{
    return MolAtomID(*this, other);
}

/** Syntactic sugar for operator+ */
MolAtomID MolID::operator&&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
MolAtomID MolID::operator&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Search for matching molecules using this ID, or other */
IDOrSet<MolID> MolID::operator*(const MolID &other) const
{
    return IDOrSet<MolID>(*this, other);
}

/** Syntactic sugar for operator* */
IDOrSet<MolID> MolID::operator||(const MolID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<MolID> MolID::operator|(const MolID &other) const
{
    return this->operator*(other);
}

/** Search for matching molecules using this ID, or other */
IDOrSet<AtomID> MolID::operator*(const AtomID &other) const
{
    return other * *this;
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> MolID::operator||(const AtomID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> MolID::operator|(const AtomID &other) const
{
    return this->operator*(other);
}

void MolID::processMatches(QList<MolNum> &matches, const Molecules &mols) const
{
    if (matches.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "There were no molecules that matched the ID \"%1\".")
                .arg(this->toString()), CODELOC );
}

///////
/////// Implementation of MolIdx
///////

static const RegisterMetaType<MolIdx> r_molidx;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MolIdx &molidx)
{
    writeHeader(ds, r_molidx, 1);

    ds << static_cast<const SireID::Index_T_<MolIdx>&>(molidx);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MolIdx &molidx)
{
    VersionID v = readHeader(ds, r_molidx);

    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<MolIdx>&>(molidx);
    }
    else
        throw version_error( v, "1", r_molidx, CODELOC );

    return ds;
}

MolIdx::MolIdx() : SireID::Index_T_<MolIdx>(), MolID()
{}

MolIdx::MolIdx(qint32 idx) : SireID::Index_T_<MolIdx>(idx), MolID()
{}

MolIdx::MolIdx(const MolIdx &other) : SireID::Index_T_<MolIdx>(other), MolID(other)
{}

MolIdx::~MolIdx()
{}

MolIdx MolIdx::null()
{
    return MolIdx();
}

bool MolIdx::isNull() const
{
    return SireID::Index_T_<MolIdx>::isNull();
}

uint MolIdx::hash() const
{
    return SireID::Index_T_<MolIdx>::hash();
}

QString MolIdx::toString() const
{
    return QString("MolIdx(%1)").arg(_idx);
}

MolIdx& MolIdx::operator=(const MolIdx &other)
{
    SireID::IndexBase::operator=(other);
    MolID::operator=(other);
    return *this;
}

bool MolIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolIdx>(*this, other);
}

QList<MolNum> MolIdx::map(const Molecules &molecules) const
{
    int i = SireID::Index(*this).map( molecules.count() );

    QList<MolNum> molnums;

    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        if (i == 0)
        {
            molnums.append(it.key());
            break;
        }
    }

    BOOST_ASSERT( not molnums.isEmpty() );

    return molnums;
}

QList<MolNum> MolIdx::map(const SelectorMol &molecules) const
{
    auto mol = molecules.molecule(this->value());
    QList<MolNum> molnums;
    molnums.append(mol.number());
    return molnums;
}

QList<MolNum> MolIdx::map(const MoleculeGroup &molgroup) const
{
    return molgroup.map(*this);
}

QList<MolNum> MolIdx::map(const MolGroupsBase &molgroups) const
{
    return molgroups.map(*this);
}

const char* MolIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolIdx>() );
}

///////
/////// Implementation of MolNum
///////

static const RegisterMetaType<MolNum> r_molnum;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MolNum &molnum)
{
    writeHeader(ds, r_molnum, 1);

    ds << static_cast<const SireID::Number&>(molnum);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MolNum &molnum)
{
    VersionID v = readHeader(ds, r_molnum);

    if (v == 1)
    {
        ds >> static_cast<SireID::Number&>(molnum);
    }
    else
        throw version_error( v, "1", r_molnum, CODELOC );

    return ds;
}

MolNum::MolNum() : SireID::Number(), MolID()
{}

MolNum::MolNum(quint32 num) : SireID::Number(num), MolID()
{}

MolNum::MolNum(const MolNum &other) : SireID::Number(other), MolID(other)
{}

MolNum::~MolNum()
{}

bool MolNum::isNull() const
{
    return SireID::Number::isNull();
}

uint MolNum::hash() const
{
    return ::qHash( static_cast<const SireID::Number&>(*this) );
}

QString MolNum::toString() const
{
    return QString("MolNum(%1)").arg(_num);
}

MolNum& MolNum::operator=(const MolNum &other)
{
    SireID::Number::operator=(other);
    MolID::operator=(other);
    return *this;
}

bool MolNum::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolNum>(*this, other);
}

bool MolNum::operator==(const MolNum &other) const
{
    return _num == other._num;
}

bool MolNum::operator!=(const MolNum &other) const
{
    return _num != other._num;
}

bool MolNum::operator<(const MolNum &other) const
{
    return _num < other._num;
}

bool MolNum::operator<=(const MolNum &other) const
{
    return _num <= other._num;
}

bool MolNum::operator>(const MolNum &other) const
{
    return _num > other._num;
}

bool MolNum::operator>=(const MolNum &other) const
{
    return _num >= other._num;
}

QList<MolNum> MolNum::map(const SelectorMol &molecules) const
{
    if (not molecules.contains(*this))
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in the set of molecules.")
                .arg(this->toString()), CODELOC );

    QList<MolNum> molnums;
    molnums.append(*this);

    return molnums;
}

QList<MolNum> MolNum::map(const Molecules &molecules) const
{
    if (not molecules.contains(*this))
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in the set of molecules.")
                .arg(this->toString()), CODELOC );

    QList<MolNum> molnums;
    molnums.append(*this);

    return molnums;
}

QList<MolNum> MolNum::map(const MoleculeGroup &molgroup) const
{
    return molgroup.map(*this);
}

QList<MolNum> MolNum::map(const MolGroupsBase &molgroups) const
{
    return molgroups.map(*this);
}

const char* MolNum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolNum>() );
}

///////
/////// Implementation of MolName
///////

static const RegisterMetaType<MolName> r_molname;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MolName &molname)
{
    writeHeader(ds, r_molname, 1);

    ds << static_cast<const SireID::Name&>(molname);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MolName &molname)
{
    VersionID v = readHeader(ds, r_molname);

    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(molname);
    }
    else
        throw version_error( v, "1", r_molname, CODELOC );

    return ds;
}

MolName::MolName() : SireID::Name(), MolID()
{}

MolName::MolName(const QString &name) : SireID::Name(name), MolID()
{}

MolName::MolName(const QString &name, SireID::CaseSensitivity case_sensitivity)
        : SireID::Name(name, case_sensitivity), MolID()
{}

MolName::MolName(const MolName &other) : SireID::Name(other), MolID(other)
{}

MolName::~MolName()
{}

bool MolName::isNull() const
{
    return SireID::Name::isNull();
}

uint MolName::hash() const
{
    return qHash(_name);
}

QString MolName::toString() const
{
    if (case_sensitive)
        return QString("MolName('%1')").arg(_name);
    else
        return QString("MolName('%1', isCaseSensitive=False)").arg(_name);
}

MolName& MolName::operator=(const MolName &other)
{
    SireID::Name::operator=(other);
    MolID::operator=(other);
    return *this;
}

bool MolName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolName>(*this, other);
}

bool MolName::operator==(const MolName &other) const
{
    return SireID::Name::operator==(other);
}

bool MolName::operator!=(const MolName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<MolNum> MolName::map(const Molecules &molecules) const
{
    QList<MolNum> molnums;

    if (this->isCaseSensitive())
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (it.value().name() == *this)
                molnums.append( it.key() );
        }
    }
    else
    {
        QString lower_name = QString(*this).toLower();

        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (QString(it.value().name()).toLower() == lower_name)
                molnums.append( it.key() );
        }
    }

    if (molnums.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with name \"%1\" in the set of molecules.")
                .arg(_name), CODELOC );

    return molnums;
}

QList<MolNum> MolName::map(const MoleculeGroup &molgroup) const
{
    return molgroup.map(*this);
}

QList<MolNum> MolName::map(const MolGroupsBase &molgroups) const
{
    return molgroups.map(*this);
}

const char* MolName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolName>() );
}

//////
////// Implementation of IDAndSet<MolID>
//////

static const RegisterMetaType< IDAndSet<MolID> > r_idandset_molid;

namespace SireID
{

/** Null constructor */
IDAndSet<MolID>::IDAndSet() : MolID()
{}

/** Add the passed ID to the list */
void IDAndSet<MolID>::add(const MolID &id)
{
    if (id.isNull())
        return;

    else if (id.isA<MolIdentifier>())
    {
        this->add(id.asA<MolIdentifier>().base());
    }
    else if (id.isA< IDAndSet<MolID> >())
        ids += id.asA< IDAndSet<MolID> >().ids;
    else
        ids.insert( MolIdentifier(id) );
}

/** Construct from the passed ID */
IDAndSet<MolID>::IDAndSet(const MolID &id) : MolID()
{
    this->add(id);
}

/** Construct from the passed IDs */
IDAndSet<MolID>::IDAndSet(const MolID &id0, const MolID &id1) : MolID()
{
    this->add(id0);
    this->add(id1);
}

/** Construct from the passed list of IDs */
IDAndSet<MolID>::IDAndSet(const QList<MolIdentifier> &new_ids) : MolID()
{
    for (QList<MolIdentifier>::const_iterator it = new_ids.constBegin();
         it != new_ids.constEnd();
         ++it)
    {
        this->add(it->base());
    }
}

/** Copy constructor */
IDAndSet<MolID>::IDAndSet(const IDAndSet &other) : MolID(other), ids(other.ids)
{}

/** Destructor */
IDAndSet<MolID>::~IDAndSet()
{}

/** Is this selection null? */
bool IDAndSet<MolID>::isNull() const
{
    return ids.isEmpty();
}

/** Return a hash of this identifier */
uint IDAndSet<MolID>::hash() const
{
    uint h = 0;

    for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        h += it->hash();
    }

    return h;
}

/** Return a string representatio of this ID */
QString IDAndSet<MolID>::toString() const
{
    if (ids.isEmpty())
        return QObject::tr("null");
    else
    {
        QStringList idstrings;

        for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
             it != ids.constEnd();
             ++it)
        {
            idstrings.append( it->toString() );
        }

        return idstrings.join( QObject::tr(" and ") );
    }
}

/** Return all of the IDs in this set */
const QSet<MolIdentifier>& IDAndSet<MolID>::IDs() const
{
    return ids;
}

/** Copy assignment operator */
IDAndSet<MolID>& IDAndSet<MolID>::operator=(const IDAndSet<MolID> &other)
{
    ids = other.ids;
    return *this;
}

/** Copy assignment operator */
IDAndSet<MolID>& IDAndSet<MolID>::operator=(const MolID &other)
{
    ids.clear();
    this->add(other);

    return *this;
}

/** Comparison operator */
bool IDAndSet<MolID>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< IDAndSet<MolID> >(*this, other);
}

/** Comparison operator */
bool IDAndSet<MolID>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

/** Comparison operator */
bool IDAndSet<MolID>::operator==(const IDAndSet<MolID> &other) const
{
    return ids == other.ids;
}

/** Comparison operator */
bool IDAndSet<MolID>::operator!=(const IDAndSet<MolID> &other) const
{
    return ids != other.ids;
}

/** Comparison operator */
bool IDAndSet<MolID>::operator==(const MolID &other) const
{
    return this->operator==( IDAndSet<MolID>(other) );
}

/** Comparison operator */
bool IDAndSet<MolID>::operator!=(const MolID &other) const
{
    return this->operator!=( IDAndSet<MolID>(other) );
}

template<class T>
QList<MolNum> IDAndSet<MolID>::_pvt_map(const T &group) const
{
    if (ids.isEmpty())
        return MolIdentifier().map(group);

    QSet<MolNum> molnums;

    QSet<MolIdentifier>::const_iterator it = ids.constBegin();

    try
    {
        molnums = convert_to_qset(it->map(group));
    }
    catch(...)
    {
        //no match
    }

    for ( ++it; it != ids.constEnd(); ++it )
    {
        if (molnums.isEmpty())
            break;

        try
        {
            molnums.intersect( convert_to_qset(it->map(group)) );
        }
        catch(...)
        {
            //no match
            molnums.clear();
        }
    }

    if (molnums.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "No molecule matches the ID \"%1\".")
                .arg(this->toString()), CODELOC );

    return molnums.values();
}

/** Map this ID to the list of indicies that match this ID

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> IDAndSet<MolID>::map(const Molecules &mols) const
{
    return this->_pvt_map(mols);
}

/** Map this ID to the list of indicies that match this ID

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> IDAndSet<MolID>::map(const MoleculeGroup &molgroup) const
{
    return this->_pvt_map(molgroup);
}

/** Map this ID to the list of indicies that match this ID

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> IDAndSet<MolID>::map(const MolGroupsBase &molgroups) const
{
    return this->_pvt_map(molgroups);
}

const char* IDAndSet<MolID>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< IDAndSet<MolID> >() );
}

} // end of namespace SireID

//////
////// Implementation of IDOrSet<MolID>
//////

static const RegisterMetaType< IDOrSet<MolID> > r_idorset_molid;

namespace SireID
{

/** Null constructor */
IDOrSet<MolID>::IDOrSet() : MolID()
{}

/** Add the passed ID to the list */
void IDOrSet<MolID>::add(const MolID &id)
{
    if (id.isNull())
        return;

    else if (id.isA<MolIdentifier>())
    {
        this->add(id.asA<MolIdentifier>().base());
    }
    else if (id.isA< IDOrSet<MolID> >())
        ids += id.asA< IDOrSet<MolID> >().ids;
    else
        ids.insert( MolIdentifier(id) );
}

/** Construct from the passed ID */
IDOrSet<MolID>::IDOrSet(const MolID &id) : MolID()
{
    this->add(id);
}

/** Construct from the passed IDs */
IDOrSet<MolID>::IDOrSet(const MolID &id0, const MolID &id1) : MolID()
{
    this->add(id0);
    this->add(id1);
}

/** Construct from the passed list of IDs */
IDOrSet<MolID>::IDOrSet(const QList<MolIdentifier> &new_ids) : MolID()
{
    for (QList<MolIdentifier>::const_iterator it = new_ids.constBegin();
         it != new_ids.constEnd();
         ++it)
    {
        this->add(it->base());
    }
}

/** Copy constructor */
IDOrSet<MolID>::IDOrSet(const IDOrSet &other) : MolID(other), ids(other.ids)
{}

/** Destructor */
IDOrSet<MolID>::~IDOrSet()
{}

/** Is this selection null? */
bool IDOrSet<MolID>::isNull() const
{
    return ids.isEmpty();
}

/** Return a hash of this identifier */
uint IDOrSet<MolID>::hash() const
{
    uint h = 0;

    for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        h += it->hash();
    }

    return h;
}

/** Return a string representatio of this ID */
QString IDOrSet<MolID>::toString() const
{
    if (ids.isEmpty())
        return QObject::tr("null");
    else
    {
        QStringList idstrings;

        for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
             it != ids.constEnd();
             ++it)
        {
            idstrings.append( it->toString() );
        }

        return idstrings.join( QObject::tr(" and ") );
    }
}

/** Return all of the IDs in this set */
const QSet<MolIdentifier>& IDOrSet<MolID>::IDs() const
{
    return ids;
}

/** Copy assignment operator */
IDOrSet<MolID>& IDOrSet<MolID>::operator=(const IDOrSet<MolID> &other)
{
    ids = other.ids;
    return *this;
}

/** Copy assignment operator */
IDOrSet<MolID>& IDOrSet<MolID>::operator=(const MolID &other)
{
    ids.clear();
    this->add(other);

    return *this;
}

/** Comparison operator */
bool IDOrSet<MolID>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< IDOrSet<MolID> >(*this, other);
}

/** Comparison operator */
bool IDOrSet<MolID>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

/** Comparison operator */
bool IDOrSet<MolID>::operator==(const IDOrSet<MolID> &other) const
{
    return ids == other.ids;
}

/** Comparison operator */
bool IDOrSet<MolID>::operator!=(const IDOrSet<MolID> &other) const
{
    return ids != other.ids;
}

/** Comparison operator */
bool IDOrSet<MolID>::operator==(const MolID &other) const
{
    return this->operator==( IDOrSet<MolID>(other) );
}

/** Comparison operator */
bool IDOrSet<MolID>::operator!=(const MolID &other) const
{
    return this->operator!=( IDOrSet<MolID>(other) );
}

QList<MolNum> IDOrSet<MolID>::process(QList<MolNum> molnums) const
{

    QSet<MolNum> set;
    set.reserve(molnums.count());

    QMutableListIterator<MolNum> it(molnums);

    while (it.hasNext())
    {
        it.next();

        if (set.contains(it.value()))
            it.remove();
        else
            set.insert(it.value());
    }

    if (molnums.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "There are no molecules that match the ID \"%1\".")
                .arg(this->toString()), CODELOC );

    return molnums;
}

/** Map this ID to the list of indicies that match this ID

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> IDOrSet<MolID>::map(const Molecules &mols) const
{
    if (ids.isEmpty())
        return MolIdentifier().map(mols);

    QList<MolNum> molnums;

    for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        try
        {
            molnums += it->map(mols);
        }
        catch(...)
        {
            //no match
        }
    }

    return this->process(molnums);
}

/** Map this ID to the list of indicies that match this ID

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> IDOrSet<MolID>::map(const MoleculeGroup &molgroup) const
{
    if (ids.isEmpty())
        return MolIdentifier().map(molgroup);

    QList<MolNum> molnums;

    for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        try
        {
            molnums += it->map(molgroup);
        }
        catch(...)
        {
            //no match
        }
    }

    return this->process(molnums);
}

/** Map this ID to the list of indicies that match this ID

    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
QList<MolNum> IDOrSet<MolID>::map(const MolGroupsBase &molgroups) const
{
    if (ids.isEmpty())
        return MolIdentifier().map(molgroups);

    QList<MolNum> molnums;

    for (QSet<MolIdentifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        try
        {
            molnums += it->map(molgroups);
        }
        catch(...)
        {
            //no match
        }
    }

    return this->process(molnums);
}

const char* IDOrSet<MolID>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< IDOrSet<MolID> >() );
}

IDOrSet<MolID>* IDOrSet<MolID>::clone() const
{
    return new IDOrSet<MolID>(*this);
}

IDAndSet<MolID>* IDAndSet<MolID>::clone() const
{
    return new IDAndSet<MolID>(*this);
}

} // end of namespace SireID

MolIdx* MolIdx::clone() const
{
    return new MolIdx(*this);
}


MolName* MolName::clone() const
{
    return new MolName(*this);
}


MolNum* MolNum::clone() const
{
    return new MolNum(*this);
}

