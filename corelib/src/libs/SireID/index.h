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

#ifndef SIREID_INDEX_H
#define SIREID_INDEX_H

#include "id.h"

#include <limits>

SIRE_BEGIN_HEADER

namespace SireID
{
class IndexBase;
class Index;
}

SIREID_EXPORT QDataStream& operator<<(QDataStream&, const SireID::IndexBase&);
SIREID_EXPORT QDataStream& operator>>(QDataStream&, SireID::IndexBase&);

SIREID_EXPORT QDataStream& operator<<(QDataStream&, const SireID::Index&);
SIREID_EXPORT QDataStream& operator>>(QDataStream&, SireID::Index&);

namespace SireID
{

/** This is the base class of all Index objects. An Index object
    provides the index of an object in an indexable list or array (or indeed
    any container that holds objects in a linear, numerical indexed
    manner (e.g. atoms in a Molecule, Molecules in a group)
    
    This class cannot be instantiated on its own - it must be 
    inherited by a derived class to be used.
    
    @author Christopher Woods
*/
class SIREID_EXPORT IndexBase
{

friend SIREID_EXPORT QDataStream& ::operator<<(QDataStream&, const IndexBase&);
friend SIREID_EXPORT QDataStream& ::operator>>(QDataStream&, IndexBase&);

public:
    ~IndexBase();

    static qint32 null();

    bool isNull() const;

    operator qint32() const;

    qint32 value() const;

    qint32 map(qint32 n) const;

    uint hash() const;

protected:
    explicit IndexBase(qint32 idx = IndexBase::null());
    
    IndexBase(const IndexBase &other);
    
    IndexBase& operator=(const IndexBase &other);
    
    void throwInvalidIndex(qint32 n) const;
    
    /** The actual index value */
    qint32 _idx;
};

/** This derived version of index provides all of the
    standard operators that you would expect.
    
    @author Christopher Woods
*/
template<class T>
class Index_T_ : public IndexBase
{
public:
    explicit Index_T_(qint32 idx=IndexBase::null());
    
    explicit Index_T_(const Index_T_<T> &other);

    ~Index_T_();

    const char* what() const;
    
    IndexBase* clone() const;
    
    uint hash() const;

    bool operator==(const T &other) const;
    bool operator==(const ID &other) const;
    bool operator==(qint32 val) const;
    bool operator!=(const T &other) const;
    bool operator!=(qint32 val) const;
    
    T& operator=(qint32 idx);
    
    T& operator+=(qint32 val);
    T& operator++();
    T operator++(qint32);
    
    T& operator-=(qint32 val);
    T& operator--();
    T operator--(qint32);
};

class SIREID_EXPORT Index : public Index_T_<Index>
{
public:
    explicit Index(qint32 idx = IndexBase::null());
    
    Index(const Index &other);
    
    ~Index();
    
    static const char* typeName();
    
    Index* clone() const;
    
    static Index null();
    
    QString toString() const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return a hash of this index */
inline uint qHash(const IndexBase &index)
{
    return index.hash();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Index_T_<T>::Index_T_(qint32 idx) : IndexBase(idx)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Index_T_<T>::Index_T_(const Index_T_<T> &other) : IndexBase(other)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Index_T_<T>::~Index_T_()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* Index_T_<T>::what() const
{
    return T::typeName();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
IndexBase* Index_T_<T>::clone() const
{
    return new T( static_cast<const T&>(*this) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
uint Index_T_<T>::hash() const
{
    return IndexBase::hash();
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Index_T_<T>::operator==(const T &other) const
{
    return _idx == other._idx;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Index_T_<T>::operator==(const ID &other) const
{
    return ID::compare<T>(static_cast<const T&>(*this),other);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Index_T_<T>::operator==(qint32 val) const
{
    return _idx == val;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Index_T_<T>::operator!=(const T &other) const
{
    return _idx != other._idx;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Index_T_<T>::operator!=(qint32 val) const
{
    return _idx != val;
}

/** Assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& Index_T_<T>::operator=(qint32 idx)
{
    _idx = idx;
    return static_cast<T&>(*this);
}

/** Increment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& Index_T_<T>::operator+=(qint32 val)
{
    _idx += val;
    return static_cast<T&>(*this);
}

/** Increment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& Index_T_<T>::operator++()
{
    ++_idx;
    return static_cast<T&>(*this);
}

/** Increment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T Index_T_<T>::operator++(qint32)
{
    T orig(*this);
    ++_idx;
    
    return orig;
}

/** Decrement operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& Index_T_<T>::operator-=(qint32 val)
{
    _idx -= val;
    return static_cast<T&>(*this);
}

/** Decrement operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& Index_T_<T>::operator--()
{
    --_idx;
    return static_cast<T&>(*this);
}

/** Decrement operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T Index_T_<T>::operator--(qint32)
{
    T orig(*this);
    --_idx;
    return orig;
}


#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireID::Index )

SIRE_EXPOSE_CLASS( SireID::IndexBase )
SIRE_EXPOSE_CLASS( SireID::Index )

SIRE_END_HEADER

#endif
