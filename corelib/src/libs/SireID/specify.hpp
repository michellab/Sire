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

#ifndef SIREID_SPECIFY_HPP
#define SIREID_SPECIFY_HPP

#include <QString>

#include "SireID/index.h"
#include "SireBase/range.h"

SIRE_BEGIN_HEADER

namespace SireID
{
template<class ID>
class Specify;
}

template<class ID>
QDataStream& operator<<(QDataStream&, const SireID::Specify<ID>&);
template<class ID>
QDataStream& operator>>(QDataStream&, SireID::Specify<ID>&);

namespace SireID
{

/** This class is used to help form specified ID matches, 
    e.g. the third residue called alanine ( ResName("ALA")[2] )
    or the last three atoms called "CA" ( AtomName("CA")(-3,-1) )
    
    @author Christopher Woods
*/
template<class ID>
class Specify : public ID
{

friend QDataStream& ::operator<<<>(QDataStream&, const Specify<ID>&);
friend QDataStream& ::operator>><>(QDataStream&, Specify<ID>&);

public:
    Specify();
    Specify(const ID &id, qint64 index);
    Specify(const ID &id, qint64 start, qint64 end);
    Specify(const ID &id, qint64 start, qint64 end, qint64 increment);
    Specify(const ID &id, const SireBase::Range &range);
    
    Specify(const Specify<ID> &other);
    
    ~Specify();
    
    static const char* typeName();
    
    const char* what() const;
    
    Specify<ID>* clone() const;

    Specify<ID>& operator=(const Specify<ID> &other);
    
    bool operator==(const Specify<ID> &other) const;
    bool operator==(const SireID::ID &other) const;

    bool operator!=(const Specify<ID> &other) const;
    
    bool operator!=(const SireID::ID &other) const;
    
    Specify<ID> operator[](qint64 i) const;
    Specify<ID> operator[](const SireBase::Range &range) const;
    
    Specify<ID> operator()(qint64 i) const;
    
    Specify<ID> operator()(qint64 start, qint64 end) const;
    Specify<ID> operator()(qint64 start, qint64 end, qint64 increment) const;
    Specify<ID> operator()(const SireBase::Range &range) const;
    
    uint hash() const;
    
    bool isNull() const;
    
    QString toString() const;
    
    QList<typename ID::Index> map(const typename ID::SearchObject &obj) const;

private:
    #ifndef SIRE_SKIP_INLINE_FUNCTIONS
    typename ID::Identifier id;
    #endif //SIRE_SKIP_INLINE_FUNCTIONS

    /** The range of indicies that are specified */
    SireBase::RangePtr idxs;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::Specify() : ID()
{}

/** Construct, using the passed ID and index */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::Specify(const ID &idobj, qint64 i)
            : ID(), id(idobj), idxs(SireBase::Range::create(i))
{}

/** Construct using the passed ID and range */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::Specify(const ID &idobj, qint64 start, qint64 end)
            : ID(), id(idobj), idxs(SireBase::Range::create(start,end))
{}
  
/** Construct using the passed ID and range */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::Specify(const ID &idobj, qint64 start, qint64 end, qint64 increment)
            : ID(), id(idobj), idxs(SireBase::Range::create(start,end,increment))
{}
  
/** Construct using the passed ID and range */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::Specify(const ID &idobj, const SireBase::Range &range)
            : ID(), id(idobj), idxs(range)
{}
  
/** Copy constructor */  
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::Specify(const Specify<ID> &other)
            : ID(other), id(other.id), idxs(other.idxs)
{}
  
/** Destructor */  
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>::~Specify()
{}

/** Copy assignment operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>& Specify<ID>::operator=(const Specify<ID> &other)
{
    if (&other != this)
    {
        ID::operator=(other);
        id = other.id;
        idxs = other.idxs;
    }
    
    return *this;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool Specify<ID>::operator==(const Specify<ID> &other) const
{
    return idxs == other.idxs and
           id == other.id;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool Specify<ID>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< Specify<ID> >(*this, other);
}

/** Return a string representation of this ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QString Specify<ID>::toString() const
{
    return QString("(%1)[%2]").arg(id.toString()).arg(idxs.read().toString());
}

/** Map this ID to the indicies of the matching objects in 'obj' */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QList<typename ID::Index> Specify<ID>::map(const typename ID::SearchObject &obj) const
{
    //first get all of the matches
    QList<typename ID::Index> found_idxs = id.map(obj);
    
    //now get the specified matches
    auto range = idxs.read().populate(found_idxs.count());

    QList<typename ID::Index> specified_idxs;
    
    while (range.read().hasNext())
    {
        specified_idxs.append( found_idxs[range.edit().next()] );
    }
    
    return specified_idxs;
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* Specify<ID>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Specify<ID> >() );
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* Specify<ID>::what() const
{
    return Specify<ID>::typeName();
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID>* Specify<ID>::clone() const
{
    return new Specify<ID>(*this);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool Specify<ID>::operator!=(const Specify<ID> &other) const
{
    return not this->operator==(other);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool Specify<ID>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

/** Return the ith item matching the ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID> Specify<ID>::operator[](qint64 i) const
{
    return Specify<ID>(*this, i);
}

/** Return the items matching the ID range */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID> Specify<ID>::operator[](const SireBase::Range &range) const
{
    return Specify<ID>(*this, range);
}

/** Return the ith item matching the ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID> Specify<ID>::operator()(qint64 i) const
{
    return Specify<ID>(*this, i);
}

/** Return the items from [start,end) matching the ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID> Specify<ID>::operator()(qint64 start, qint64 end) const
{
    return Specify<ID>(*this, start, end);
}

/** Return the items from [start,end,increment) matching the ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID> Specify<ID>::operator()(qint64 start, qint64 end, qint64 increment) const
{
    return Specify<ID>(*this, start, end, increment);
}

/** Return the items matching the ID range */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
Specify<ID> Specify<ID>::operator()(const SireBase::Range &range) const
{
    return Specify<ID>(*this, range);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
uint Specify<ID>::hash() const
{
    return id.hash() + 5;
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool Specify<ID>::isNull() const
{
    return id.isNull();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise to a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireID::Specify<ID> &id)
{
    ds << id.id << id.idxs;
    return ds;
}

/** Extract from a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireID::Specify<ID> &id)
{
    ds >> id.id >> id.idxs;
    return ds;
}

SIRE_END_HEADER

#endif
