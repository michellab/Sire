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

#include "residentifier.h"
#include "residx.h"
#include "resnum.h"
#include "resname.h"
#include "molinfo.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

////////
//////// Implementation of ResIdentifier
////////

static const RegisterMetaType<ResIdentifier> r_resid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ResIdentifier &resid)
{
    writeHeader(ds, r_resid, 1);
    
    SireStream::savePolyPointer(ds, resid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResIdentifier &resid)
{
    VersionID v = readHeader(ds, r_resid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, resid.d);
    }
    else
        throw version_error( v, "1", r_resid, CODELOC );
        
    return ds;
}

/** Null constructor */
ResIdentifier::ResIdentifier() : ResID()
{}

/** Construct from the passed ResID */
ResIdentifier::ResIdentifier(const ResID &resid)
              : ResID()
{
    if (resid.isA<ResIdentifier>())
        d = resid.asA<ResIdentifier>().d;
    else if (not resid.isNull())
        d.reset( resid.clone() );
}

/** Copy constructor */
ResIdentifier::ResIdentifier(const ResIdentifier &other)
              : ResID(other), d(other.d)
{}

/** Destructor */
ResIdentifier::~ResIdentifier()
{}

/** Is this selection null? */
bool ResIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint ResIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString ResIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const ResID& ResIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
ResIdentifier& ResIdentifier::operator=(const ResIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
ResIdentifier& ResIdentifier::operator=(const ResID &other)
{
    if (other.isA<ResIdentifier>())
        d = other.asA<ResIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool ResIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ResIdentifier>(*this, other);
}

/** Comparison operator */
bool ResIdentifier::operator==(const ResIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool ResIdentifier::operator!=(const ResIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool ResIdentifier::operator==(const ResID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<ResIdentifier>())
        return this->operator==(other.asA<ResIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool ResIdentifier::operator!=(const ResID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<ResIdentifier>())
        return this->operator!=(other.asA<ResIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID to the list of indicies of reidues that match this ID

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> ResIdentifier::map(const MolInfo &molinfo) const
{
    if (d.get() == 0)
        return molinfo.getResidues();
    else
        return d->map(molinfo);
}

const char* ResIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResIdentifier>() );
}

///////
/////// Implementation of ResIdx
///////

static const RegisterMetaType<ResIdx> r_residx;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ResIdx &residx)
{
    writeHeader(ds, r_residx, 1);
    
    ds << static_cast<const SireID::Index_T_<ResIdx>&>(residx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResIdx &residx)
{
    VersionID v = readHeader(ds, r_residx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<ResIdx>&>(residx);
    }
    else
        throw version_error( v, "1", r_residx, CODELOC );
        
    return ds;
}

ResIdx::ResIdx() : SireID::Index_T_<ResIdx>(), ResID()
{}

ResIdx::ResIdx(quint32 idx) 
          : SireID::Index_T_<ResIdx>(idx), ResID()
{}

ResIdx::ResIdx(const ResIdx &other) 
          : SireID::Index_T_<ResIdx>(other), ResID(other)
{}

ResIdx::~ResIdx()
{}

ResIdx ResIdx::null()
{
    return ResIdx();
}

bool ResIdx::isNull() const
{
    return SireID::Index_T_<ResIdx>::isNull();
}

uint ResIdx::hash() const
{
    return SireID::Index_T_<ResIdx>::hash();
}

QString ResIdx::toString() const
{
    return QString("ResIdx(%1)").arg(_idx);
}

ResIdx& ResIdx::operator=(const ResIdx &other)
{
    SireID::IndexBase::operator=(other);
    ResID::operator=(other);
    return *this;
}

bool ResIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ResIdx>(*this, other);
}

QList<ResIdx> ResIdx::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* ResIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResIdx>() );
}

///////
/////// Implementation of ResName
///////

static const RegisterMetaType<ResName> r_resname;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ResName &resname)
{
    writeHeader(ds, r_resname, 1);
    
    ds << static_cast<const SireID::Name&>(resname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResName &resname)
{
    VersionID v = readHeader(ds, r_resname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(resname);
    }
    else
        throw version_error( v, "1", r_resname, CODELOC );
        
    return ds;
}

ResName::ResName() : SireID::Name(), ResID()
{}

ResName::ResName(const QString &name) : SireID::Name(name), ResID()
{}

ResName::ResName(const QString &name, SireID::CaseSensitivity case_sensitivity)
        : SireID::Name(name, case_sensitivity), ResID()
{}

ResName::ResName(const ResName &other) : SireID::Name(other), ResID(other)
{}

ResName::~ResName()
{}

bool ResName::isNull() const
{
    return SireID::Name::isNull();
}

uint ResName::hash() const
{
    return qHash(_name);
}

QString ResName::toString() const
{
    if (case_sensitive)
        return QString("ResName('%1')").arg(_name);
    else
        return QString("ResName('%1', isCaseSensitive=False)").arg(_name);
}

ResName& ResName::operator=(const ResName &other)
{
    SireID::Name::operator=(other);
    ResID::operator=(other);
    return *this;
}

bool ResName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ResName>(*this, other);
}

bool ResName::operator==(const ResName &other) const
{
    return SireID::Name::operator==(other);
}

bool ResName::operator!=(const ResName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<ResIdx> ResName::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* ResName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResName>() );
}

///////
/////// Implementation of ResNum
///////

static const RegisterMetaType<ResNum> r_resnum;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ResNum &resnum)
{
    writeHeader(ds, r_resnum, 1);
    
    ds << static_cast<const SireID::Number&>(resnum);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResNum &resnum)
{
    VersionID v = readHeader(ds, r_resnum);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Number&>(resnum);
    }
    else
        throw version_error( v, "1", r_resnum, CODELOC );
        
    return ds;
}

ResNum::ResNum() : SireID::Number(), ResID()
{}

ResNum::ResNum(quint32 num) : SireID::Number(num), ResID()
{}

ResNum::ResNum(const ResNum &other) : SireID::Number(other), ResID(other)
{}

ResNum::~ResNum()
{}

bool ResNum::isNull() const
{
    return SireID::Number::isNull();
}

uint ResNum::hash() const
{
    return ::qHash( static_cast<const SireID::Number&>(*this) );
}

QString ResNum::toString() const
{
    return QString("ResNum(%1)").arg(_num);
}

ResNum& ResNum::operator=(const ResNum &other)
{
    SireID::Number::operator=(other);
    ResID::operator=(other);
    return *this;
}

bool ResNum::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ResNum>(*this, other);
}

bool ResNum::operator==(const ResNum &other) const
{
    return _num == other._num;
}

bool ResNum::operator!=(const ResNum &other) const
{
    return _num != other._num;
}

bool ResNum::operator<(const ResNum &other) const
{
    return _num < other._num;
}

bool ResNum::operator<=(const ResNum &other) const
{
    return _num <= other._num;
}

bool ResNum::operator>(const ResNum &other) const
{
    return _num > other._num;
}

bool ResNum::operator>=(const ResNum &other) const
{
    return _num >= other._num;
}

QList<ResIdx> ResNum::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* ResNum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResNum>() );
}

ResIdentifier* ResIdentifier::clone() const
{
    return new ResIdentifier(*this);
}

ResName* ResName::clone() const
{
    return new ResName(*this);
}


ResNum* ResNum::clone() const
{
    return new ResNum(*this);
}


ResIdx* ResIdx::clone() const
{
    return new ResIdx(*this);
}

