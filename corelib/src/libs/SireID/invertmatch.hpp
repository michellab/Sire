/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREID_INVERTMATCH_HPP
#define SIREID_INVERTMATCH_HPP

#include <QString>

#include "SireID/index.h"

SIRE_BEGIN_HEADER

namespace SireID
{
template<class ID>
class InvertMatch;
}

template<class ID>
QDataStream& operator<<(QDataStream&, const SireID::InvertMatch<ID>&);
template<class ID>
QDataStream& operator>>(QDataStream&, SireID::InvertMatch<ID>&);

namespace SireID
{

/** This class is used to form the inverse or "not match"
    of any ID, e.g. ResName("GLY").inverse() will match
    any residue that is not GLY.
    
    @author Christopher Woods
*/
template<class ID>
class InvertMatch : public ID
{

friend QDataStream& ::operator<<<>(QDataStream&, const InvertMatch<ID>&);
friend QDataStream& ::operator>><>(QDataStream&, InvertMatch<ID>&);

public:
    InvertMatch();
    InvertMatch(const ID &id);
    
    InvertMatch(const InvertMatch<ID> &other);
    
    ~InvertMatch();
    
    static const char* typeName();
    
    const char* what() const;
    
    InvertMatch<ID>* clone() const;

    InvertMatch<ID>& operator=(const InvertMatch<ID> &other);
    
    bool operator==(const InvertMatch<ID> &other) const;
    bool operator==(const SireID::ID &other) const;

    bool operator!=(const InvertMatch<ID> &other) const;
    
    bool operator!=(const SireID::ID &other) const;
    
    uint hash() const;
    
    bool isNull() const;
    
    QString toString() const;
    
    QList<typename ID::Index> map(const typename ID::SearchObject &obj) const;

private:
    #ifndef SIRE_SKIP_INLINE_FUNCTIONS
    typename ID::Identifier id;
    #endif //SIRE_SKIP_INLINE_FUNCTIONS
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
InvertMatch<ID>::InvertMatch() : ID()
{}

/** Construct, using the passed ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
InvertMatch<ID>::InvertMatch(const ID &idobj)
            : ID(), id(idobj)
{}
  
/** Copy constructor */  
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
InvertMatch<ID>::InvertMatch(const InvertMatch<ID> &other)
            : ID(other), id(other.id)
{}
  
/** Destructor */  
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
InvertMatch<ID>::~InvertMatch()
{}

/** Copy assignment operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
InvertMatch<ID>& InvertMatch<ID>::operator=(const InvertMatch<ID> &other)
{
    if (&other != this)
    {
        ID::operator=(other);
        id = other.id;
    }
    
    return *this;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool InvertMatch<ID>::operator==(const InvertMatch<ID> &other) const
{
    return id == other.id;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool InvertMatch<ID>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< InvertMatch<ID> >(*this, other);
}

/** Return a string representation of this ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QString InvertMatch<ID>::toString() const
{
    return QObject::tr("not { %1 }").arg(id.toString());
}

/** Map this ID to the indicies of the matching objects in 'obj' */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QList<typename ID::Index> InvertMatch<ID>::map(const typename ID::SearchObject &obj) const
{
    //first, get all of the possible matches
    QList<typename ID::Index> all_idxs = ID::matchAll(obj);

    //now get all of the matches
    QList<typename ID::Index> idxs;
    
    try
    {
        //try to map, as if there is no match, then this
        //should match everything
        idxs = id.map(obj);
    }
    catch(...)
    {}
    
    if (not idxs.isEmpty())
    {
        //now remove any match in all_idxs that is in 'idxs'
        QMutableListIterator<typename ID::Index> it(all_idxs);
        
        while (it.hasNext())
        {
            typename ID::Index idx = it.next();
        
            if (idxs.contains(idx))
            {
                it.remove();
            }
        }
    }
    
    this->processMatches(all_idxs, obj);
    
    return all_idxs;
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* InvertMatch<ID>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< InvertMatch<ID> >() );
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* InvertMatch<ID>::what() const
{
    return InvertMatch<ID>::typeName();
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
InvertMatch<ID>* InvertMatch<ID>::clone() const
{
    return new InvertMatch<ID>(*this);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool InvertMatch<ID>::operator!=(const InvertMatch<ID> &other) const
{
    return not this->operator==(other);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool InvertMatch<ID>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
uint InvertMatch<ID>::hash() const
{
    return id.hash() + 1;
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool InvertMatch<ID>::isNull() const
{
    return id.isNull();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise to a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireID::InvertMatch<ID> &id)
{
    ds << id.id;
    return ds;
}

/** Extract from a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireID::InvertMatch<ID> &id)
{
    ds >> id.id;
    return ds;
}

SIRE_END_HEADER

#endif
