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

#include "index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireID;
using namespace SireStream;

////////
//////// Implementation of IndexBase
////////

/** Serialise an Idx class */
QDataStream SIREID_EXPORT &operator<<(QDataStream &ds, const IndexBase &idx)
{
    ds << idx._idx;
    return ds;
}

/** Deserialise an Idx class */
QDataStream SIREID_EXPORT &operator>>(QDataStream &ds, IndexBase &idx)
{
    ds >> idx._idx;
    return ds;
}

IndexBase::IndexBase(qint32 idx) : _idx(idx)
{}

IndexBase::IndexBase(const IndexBase &other) : _idx(other._idx)
{}

IndexBase& IndexBase::operator=(const IndexBase &other)
{
    _idx = other._idx;
    return *this;
}

IndexBase::~IndexBase()
{}

/** Return the null index */
qint32 IndexBase::null()
{
    return std::numeric_limits<qint32>::min();
}

/** Return whether this is a null index - a null
    index is one that equals std::numeric_limits<qint32>::min(),
    which should be -21474833648 for a 32bit integer */
bool IndexBase::isNull() const
{
    return _idx == IndexBase::null();
}

/** Hash this Index */
uint IndexBase::hash() const
{
    return quint32(_idx);
}

/** Return the raw value of this index */
qint32 IndexBase::value() const
{
    return _idx;
}

/** Allow implicit conversion back to an int */
IndexBase::operator qint32() const
{
    return _idx;
}
    
void IndexBase::throwInvalidIndex(qint32 n) const
{
    if (n == 0)
        throw SireError::invalid_index(QObject::tr(
            "Cannot access item at index %1 as the container is empty!")
                .arg(_idx), CODELOC );

    else if (n == 1)
        throw SireError::invalid_index(QObject::tr(
            "Cannot access item at index %1 as there is only item in the "
            "container.").arg(_idx), CODELOC );
    
    else
        throw SireError::invalid_index( QObject::tr(
            "No item at index %1. Index range is from %2 to %3.")
                .arg(_idx).arg(-n).arg(n-1), CODELOC );
}

/** Map this index into the container of 'n' elements - this
    maps the index (with negative indexing, e.g. -1 is the last
    element), and throws an exception if the index is out 
    of the bounds of the array
    
    \throw SireError::invalid_index
*/
qint32 IndexBase::map(qint32 n) const
{
    if (_idx >= 0 and _idx < n)
        return _idx;
    else if (_idx < 0 and _idx >= -n)
        return n + _idx;
    else
    {
        throwInvalidIndex(n);
        return null();
    }
}

////////
//////// Implementation of Index
////////

static const RegisterMetaType<Index> r_index(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREID_EXPORT &operator<<(QDataStream &ds, const Index &index)
{
    writeHeader(ds, r_index, 1);
    
    ds << static_cast<const Index_T_<Index>&>(index);
    
    return ds;
}
    
/** Extract from a binary datastream */
QDataStream SIREID_EXPORT &operator>>(QDataStream &ds, Index &index)
{
    VersionID v = readHeader(ds, r_index);
    
    if (v == 1)
    {
        ds >> static_cast<Index_T_<Index>&>(index);
    }
    else
        throw version_error( v, "1", r_index, CODELOC );
    
    return ds;
}

Index::Index(qint32 idx) : Index_T_<Index>(idx)
{}

Index::Index(const Index &other) : Index_T_<Index>(other)
{}

Index::~Index()
{}

Index Index::null()
{
    return Index();
}

QString Index::toString() const
{
    return QString("Index(%1)").arg(_idx);
}

const char* Index::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Index>() );
}
