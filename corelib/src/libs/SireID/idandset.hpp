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

#ifndef SIREID_IDANDSET_HPP
#define SIREID_IDANDSET_HPP

#include <QSet>
#include <QString>
#include <QStringList>
#include <QObject>

#include "SireID/id.h"

SIRE_BEGIN_HEADER

namespace SireID
{
template<class ID>
class IDAndSet;
}

template<class ID>
QDataStream& operator<<(QDataStream&, const SireID::IDAndSet<ID>&);
template<class ID>
QDataStream& operator>>(QDataStream&, SireID::IDAndSet<ID>&);

namespace SireID
{

/** This class holds a set of IDs, thereby allowing for
    "and" matching of IDs 
    
    @author Christopher Woods
*/
template<class ID>
class SIREID_EXPORT IDAndSet : public ID
{

friend QDataStream& ::operator<<<>(QDataStream&, const IDAndSet<ID>&);
friend QDataStream& ::operator>><>(QDataStream&, IDAndSet<ID>&);

public:
    typedef typename ID::Identifier Identifier;
    typedef typename ID::Index Index;
    typedef typename ID::SearchObject SearchObject;

    IDAndSet();
    IDAndSet(const ID &id);
    IDAndSet(const ID &id0, const ID &id1);
    
    IDAndSet(const QList<typename ID::Identifier> &ids);
    
    template<class T>
    IDAndSet(const T &ids);
    
    IDAndSet(const IDAndSet<ID> &other);
    
    ~IDAndSet();
    
    static const char* typeName();
    
    const char* what() const;
    
    IDAndSet<ID>* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;

    const QSet<Identifier>& IDs() const;
    
    IDAndSet<ID>& operator=(const IDAndSet<ID> &other);
    IDAndSet<ID>& operator=(const ID &other);
    
    bool operator==(const SireID::ID &other) const;
    bool operator!=(const SireID::ID &other) const;
   
    bool operator==(const IDAndSet<ID> &other) const;
    bool operator!=(const IDAndSet<ID> &other) const;
    
    bool operator==(const ID &other) const;
    bool operator!=(const ID &other) const;
    
    QList<Index> map(const SearchObject &obj) const;

private:
    void add(const ID &id);

    QSet<typename ID::Identifier> ids;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::IDAndSet() : ID()
{}

/** Add the passed ID to the list */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
void IDAndSet<ID>::add(const ID &id)
{
    if (id.isNull())
        return;

    else if (id.template isA<Identifier>())
    {
        this->add(id.template asA<Identifier>().base());
    }
    else if (id.template isA< IDAndSet<ID> >())
        ids += id.template asA< IDAndSet<ID> >().ids;
    else
        ids.insert( Identifier(id) );
}

/** Construct from the passed ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::IDAndSet(const ID &id) : ID()
{
    this->add(id);
}

/** Construct from the passed IDs */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::IDAndSet(const ID &id0, const ID &id1) : ID()
{
    this->add(id0);
    this->add(id1);
}

/** Construct from the passed list of IDs */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::IDAndSet(const QList<typename ID::Identifier> &new_ids) : ID()
{
    for (typename QList<Identifier>::const_iterator it = new_ids.constBegin();
         it != new_ids.constEnd();
         ++it)
    {
        this->add(it->base());
    }
}

/** Construct from the passed list of IDs */
template<class ID>
template<class T>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::IDAndSet(const T &new_ids) : ID()
{
    for (typename T::const_iterator it = new_ids.constBegin();
         it != new_ids.constEnd();
         ++it)
    {
        this->add(*it);
    }
}

/** Copy constructor */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::IDAndSet(const IDAndSet &other) : ID(other), ids(other.ids)
{}

/** Destructor */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>::~IDAndSet()
{}
    
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* IDAndSet<ID>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< IDAndSet<ID> >() );
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* IDAndSet<ID>::what() const
{
    return IDAndSet<ID>::typeName();
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>* IDAndSet<ID>::clone() const
{
    return new IDAndSet<ID>(*this);
}

/** Is this selection null? */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::isNull() const
{
    return ids.isEmpty();
}

/** Return a hash of this identifier */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
uint IDAndSet<ID>::hash() const
{
    uint h = 0;
    
    for (typename QSet<Identifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        h += it->hash();
    }
    
    return h;
}
            
/** Return a string representatio of this ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QString IDAndSet<ID>::toString() const
{
    if (ids.isEmpty())
        return QObject::tr("null");
    else
    {
        QStringList idstrings;
        
        for (typename QSet<Identifier>::const_iterator it = ids.constBegin();
             it != ids.constEnd();
             ++it)
        {
            idstrings.append( it->toString() );
        }
    
        return idstrings.join( QObject::tr(" and ") );
    }
}

/** Return all of the IDs in this set */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const QSet<typename ID::Identifier>& IDAndSet<ID>::IDs() const
{
    return ids;
}

/** Copy assignment operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>& IDAndSet<ID>::operator=(const IDAndSet<ID> &other)
{
    ids = other.ids;
    return *this;
}

/** Copy assignment operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
IDAndSet<ID>& IDAndSet<ID>::operator=(const ID &other)
{
    ids.clear();
    this->add(other);
    
    return *this;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< IDAndSet<ID> >(*this, other);
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::operator==(const IDAndSet<ID> &other) const
{
    return ids == other.ids;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::operator!=(const IDAndSet<ID> &other) const
{
    return ids != other.ids;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::operator==(const ID &other) const
{
    return this->operator==( IDAndSet<ID>(other) );
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool IDAndSet<ID>::operator!=(const ID &other) const
{
    return this->operator!=( IDAndSet<ID>(other) );
}

/** Map this ID to the list of indicies that match this ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QList<typename ID::Index> IDAndSet<ID>::map(const typename ID::SearchObject &obj) const
{
    if (ids.isEmpty())
        return Identifier().map(obj);

    typename QSet<Identifier>::const_iterator it = ids.constBegin();

    QSet<Index> idxs;
    
    try
    {
        idxs = it->map(obj).toSet();
    }
    catch(...)
    {
        //no match
    }
        
    for ( ++it; it != ids.constEnd(); ++it )
    {
        if (idxs.isEmpty())
            break;
    
        try
        {
            idxs.intersect( it->map(obj).toSet() );
        }
        catch(...)
        {
            //no match
            idxs.clear();
        }
    }
    
    QList<Index> matches = idxs.toList();
    ID::processMatches(matches, obj);

    return matches;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS
/** Serialise to a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireID::IDAndSet<ID> &idset)
{
    ds << idset.ids;
    return ds;
}

/** Extract from a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireID::IDAndSet<ID> &idset)
{
    ds >> idset.ids;
    return ds;
}
#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
