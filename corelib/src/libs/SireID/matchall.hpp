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

#ifndef SIREID_MATCHALL_HPP
#define SIREID_MATCHALL_HPP

#include <QString>

#include "SireID/index.h"

SIRE_BEGIN_HEADER

namespace SireID
{
template<class ID>
class MatchAll;
}

template<class ID>
QDataStream& operator<<(QDataStream&, const SireID::MatchAll<ID>&);
template<class ID>
QDataStream& operator>>(QDataStream&, SireID::MatchAll<ID>&);

namespace SireID
{

/** This class is used as a convenience class to held provide an
    ID that matches ALL items. This is useful when combining IDs
    of different types, e.g. to turn a ResID into an ID that
    matches all atoms in that residues you could use 
    ResName("ALA") + AtomID.any()
 
    @author Christopher Woods
*/
template<class ID>
class MatchAll : public ID
{

friend QDataStream& ::operator<<<>(QDataStream&, const MatchAll<ID>&);
friend QDataStream& ::operator>><>(QDataStream&, MatchAll<ID>&);

public:
    MatchAll();
    
    MatchAll(const MatchAll<ID> &other);
    
    ~MatchAll();
    
    static const char* typeName();
    
    const char* what() const;
    
    MatchAll<ID>* clone() const;

    MatchAll<ID>& operator=(const MatchAll<ID> &other);
    
    bool operator==(const MatchAll<ID> &other) const;
    bool operator==(const SireID::ID &other) const;

    bool operator!=(const MatchAll<ID> &other) const;
    
    bool operator!=(const SireID::ID &other) const;
    
    uint hash() const;
    
    bool isNull() const;
    
    QString toString() const;
    
    QList<typename ID::Index> map(const typename ID::SearchObject &obj) const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
MatchAll<ID>::MatchAll() : ID()
{}
  
/** Copy constructor */  
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
MatchAll<ID>::MatchAll(const MatchAll<ID> &other)
            : ID(other)
{}
  
/** Destructor */  
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
MatchAll<ID>::~MatchAll()
{}

/** Copy assignment operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
MatchAll<ID>& MatchAll<ID>::operator=(const MatchAll<ID> &other)
{
    if (&other != this)
    {
        ID::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool MatchAll<ID>::operator==(const MatchAll<ID> &other) const
{
    return true;
}

/** Comparison operator */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool MatchAll<ID>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< MatchAll<ID> >(*this, other);
}

/** Return a string representation of this ID */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QString MatchAll<ID>::toString() const
{
    return QString("%1::any").arg( QString(ID::typeName()).split("::").last() );
}

/** Map this ID to the indicies of the matching objects in 'obj' */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QList<typename ID::Index> MatchAll<ID>::map(const typename ID::SearchObject &obj) const
{
    QList<typename ID::Index> idxs = ID::matchAll(obj);
    this->processMatches(idxs, obj);
    return idxs;
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* MatchAll<ID>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< MatchAll<ID> >() );
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
const char* MatchAll<ID>::what() const
{
    return MatchAll<ID>::typeName();
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
MatchAll<ID>* MatchAll<ID>::clone() const
{
    return new MatchAll<ID>(*this);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool MatchAll<ID>::operator!=(const MatchAll<ID> &other) const
{
    return not this->operator==(other);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool MatchAll<ID>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
uint MatchAll<ID>::hash() const
{
    return 1;
}

template<class ID>
SIRE_OUTOFLINE_TEMPLATE
bool MatchAll<ID>::isNull() const
{
    return false;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise to a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireID::MatchAll<ID> &id)
{
    return ds;
}

/** Extract from a binary datastream */
template<class ID>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireID::MatchAll<ID> &id)
{
    return ds;
}

SIRE_END_HEADER

#endif
