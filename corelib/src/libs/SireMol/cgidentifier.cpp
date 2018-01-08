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

#include "cgidentifier.h"
#include "cgidx.h"
#include "cgname.h"
#include "molinfo.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

/////////
///////// Implementation of CGIdentifier
/////////

static const RegisterMetaType<CGIdentifier> r_cgid;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const CGIdentifier &cgid)
{
    writeHeader(ds, r_cgid, 1);
    
    SireStream::savePolyPointer(ds, cgid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, CGIdentifier &cgid)
{
    VersionID v = readHeader(ds, r_cgid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, cgid.d);
    }
    else
        throw version_error( v, "1", r_cgid, CODELOC );
        
    return ds;
}

/** Null constructor */
CGIdentifier::CGIdentifier() : CGID()
{}

/** Construct from the passed CGID */
CGIdentifier::CGIdentifier(const CGID &cgid)
             : CGID()
{
    if (cgid.isA<CGIdentifier>())
        d = cgid.asA<CGIdentifier>().d;
    else if (not cgid.isNull())
        d.reset( cgid.clone() );
}

/** Copy constructor */
CGIdentifier::CGIdentifier(const CGIdentifier &other)
             : CGID(other), d(other.d)
{}

/** Destructor */
CGIdentifier::~CGIdentifier()
{}

/** Is this selection null? */
bool CGIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint CGIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString CGIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const CGID& CGIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
CGIdentifier& CGIdentifier::operator=(const CGIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
CGIdentifier& CGIdentifier::operator=(const CGID &other)
{
    if (other.isA<CGIdentifier>())
        d = other.asA<CGIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool CGIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<CGIdentifier>(*this, other);
}

/** Comparison operator */
bool CGIdentifier::operator==(const CGIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool CGIdentifier::operator!=(const CGIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool CGIdentifier::operator==(const CGID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<CGIdentifier>())
        return this->operator==(other.asA<CGIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool CGIdentifier::operator!=(const CGID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<CGIdentifier>())
        return this->operator!=(other.asA<CGIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID to the list of indicies of CutGroups that match this ID

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QList<CGIdx> CGIdentifier::map(const MolInfo &molinfo) const
{
    if (d.get() == 0)
        return molinfo.getCutGroups();
    else
        return d->map(molinfo);
}

const char* CGIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGIdentifier>() );
}

/////////
///////// Implementation of CGIdx
/////////

static const RegisterMetaType<CGIdx> r_cgidx;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const CGIdx &cgidx)
{
    writeHeader(ds, r_cgidx, 1);
    
    ds << static_cast<const SireID::Index_T_<CGIdx>&>(cgidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, CGIdx &cgidx)
{
    VersionID v = readHeader(ds, r_cgidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<CGIdx>&>(cgidx);
    }
    else
        throw version_error( v, "1", r_cgidx, CODELOC );
        
    return ds;
}

CGIdx::CGIdx() : SireID::Index_T_<CGIdx>(), CGID()
{}

CGIdx::CGIdx(qint32 idx) 
          : SireID::Index_T_<CGIdx>(idx), CGID()
{}

CGIdx::CGIdx(const CGIdx &other) 
          : SireID::Index_T_<CGIdx>(other), CGID(other)
{}

CGIdx::~CGIdx()
{}

CGIdx CGIdx::null()
{
    return CGIdx();
}

bool CGIdx::isNull() const
{
    return SireID::Index_T_<CGIdx>::isNull();
}

uint CGIdx::hash() const
{
    return SireID::Index_T_<CGIdx>::hash();
}

QString CGIdx::toString() const
{
    return QString("CGIdx(%1)").arg(_idx);
}

CGIdx& CGIdx::operator=(const CGIdx &other)
{
    SireID::IndexBase::operator=(other);
    CGID::operator=(other);
    return *this;
}

bool CGIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<CGIdx>(*this, other);
}

QList<CGIdx> CGIdx::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* CGIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGIdx>() );
}

/////////
///////// Implementation of CGName
/////////

static const RegisterMetaType<CGName> r_cgname;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const CGName &cgname)
{
    writeHeader(ds, r_cgname, 1);
    
    ds << static_cast<const SireID::Name&>(cgname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, CGName &cgname)
{
    VersionID v = readHeader(ds, r_cgname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(cgname);
    }
    else
        throw version_error( v, "1", r_cgname, CODELOC );
        
    return ds;
}

CGName::CGName() : SireID::Name(), CGID()
{}

CGName::CGName(const QString &name) : SireID::Name(name), CGID()
{}

CGName::CGName(const QString &name, SireID::CaseSensitivity case_sensitivity)
       : SireID::Name(name, case_sensitivity), CGID()
{}

CGName::CGName(const CGName &other) : SireID::Name(other), CGID(other)
{}

CGName::~CGName()
{}

bool CGName::isNull() const
{
    return SireID::Name::isNull();
}

uint CGName::hash() const
{
    return ::qHash(_name);
}

QString CGName::toString() const
{
    if (case_sensitive)
        return QString("CGName('%1')").arg(_name);
    else
        return QString("CGName('%1', isCaseSensitive=False)").arg(_name);
}

CGName& CGName::operator=(const CGName &other)
{
    SireID::Name::operator=(other);
    CGID::operator=(other);
    return *this;
}

bool CGName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<CGName>(*this, other);
}

bool CGName::operator==(const CGName &other) const
{
    return SireID::Name::operator==(other);
}

bool CGName::operator!=(const CGName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<CGIdx> CGName::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* CGName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGName>() );
}

CGIdentifier* CGIdentifier::clone() const
{
    return new CGIdentifier(*this);
}


CGName* CGName::clone() const
{
    return new CGName(*this);
}


CGIdx* CGIdx::clone() const
{
    return new CGIdx(*this);
}

