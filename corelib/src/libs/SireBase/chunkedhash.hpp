/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREBASE_CHUNKEDHASH_HPP
#define SIREBASE_CHUNKEDHASH_HPP

#include "sireglobal.h"

#include <QHash>

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class Key, class T, int N>
class ChunkedHash;
}

template<class Key, class T, int N>
QDataStream& operator<<(QDataStream&, const SireBase::ChunkedHash<Key,T,N>&);
template<class Key, class T, int N>
QDataStream& operator>>(QDataStream&, SireBase::ChunkedHash<Key,T,N>&);

namespace SireBase
{

namespace detail
{
template<class Key, class T, int N>
const void* get_shared_container_pointer(const SireBase::ChunkedHash<Key,T,N>&);
} // end of namespace detail

/** This is a hash that stores the values and keys in separate
    chunks - this prevents large copies if only small parts
    of the hash are changed at a time
    
    @author Christopher Woods
*/
template<class Key, class T, int N=100>
class ChunkedHash
{

friend QDataStream& ::operator<<<>(QDataStream&, const ChunkedHash<Key,T,N>&);
friend QDataStream& ::operator>><>(QDataStream&, ChunkedHash<Key,T,N>&);

friend const void* 
SireBase::detail::get_shared_container_pointer<>(const ChunkedHash<Key,T,N>&);

public:
    class iterator;
    class const_iterator;

    /** Iterator over a ChunkedHash that is allowed to modify its contents */
    class iterator
    {
    
    friend class ChunkedHash;
    friend class const_iterator;
    
    public:
        iterator();
        iterator(const iterator &other);
        
        ~iterator();
        
        iterator& operator=(const iterator &other);
        iterator& operator=(const const_iterator &other);
        
        bool operator==(const iterator &other) const;
        bool operator!=(const iterator &other) const;
        bool operator==(const const_iterator &other) const;
        bool operator!=(const const_iterator &other) const;

        const Key& key() const;
        
        T& value() const;
        
        T& operator*() const;
        T* operator->() const;
        
        iterator operator+(int j) const;
        iterator& operator++();
        iterator operator++(int);
        iterator& operator+=(int j);
    
        iterator operator-(int j) const;
        iterator& operator--();
        iterator operator--( int);
        iterator& operator-=( int j );   

    private:
        /** Pointer to the parent's chunks */
        QVector< QHash<Key,T> > *chunks;
        
        /** Index of the current chunk */
        int current_chunk;
        
        /** Iterator over the values in the current chunk */
        typename QHash<Key,T>::iterator current_it;
    };

    /** Iterator over a ChunkedHash that is not allowed to modify its contents */
    class const_iterator
    {
    
    friend class ChunkedHash;
    friend class iterator;
    
    public:
        const_iterator();
        const_iterator(const iterator &other);
        const_iterator(const const_iterator &other);
        
        ~const_iterator();
        
        const_iterator& operator=(const iterator &other);
        const_iterator& operator=(const const_iterator &other);
        
        bool operator==(const const_iterator &other) const;
        bool operator!=(const const_iterator &other) const;
        bool operator==(const iterator &other) const;
        bool operator!=(const iterator &other) const;
        
        const Key& key() const;
        const T& value() const;
        
        const T& operator*() const;
        const T* operator->() const;
        
        const_iterator operator+(int j) const;
        const_iterator& operator++();
        const_iterator operator++(int);
        const_iterator& operator+=(int j);
        
        const_iterator operator-(int j) const;
        const_iterator& operator--();
        const_iterator operator--(int);
        const_iterator& operator-=(int j);

    private:
        /** Pointer to the parent's chunks */
        const QVector< QHash<Key,T> > *chunks;
        
        /** Index of the current chunk */
        int current_chunk;
        
        /** Iterator over the values in the current chunk */
        typename QHash<Key,T>::const_iterator current_it;
    };

    ChunkedHash();
    ChunkedHash(const QHash<Key,T> &hash);
    
    ChunkedHash(const ChunkedHash<Key,T,N> &other);
    
    ~ChunkedHash();
    
    ChunkedHash<Key,T,N>& operator=(const QHash<Key,T> &hash);
    ChunkedHash<Key,T,N>& operator=(const ChunkedHash<Key,T,N> &other);
    
    bool operator==(const ChunkedHash<Key,T,N> &other) const;
    bool operator!=(const ChunkedHash<Key,T,N> &other) const;
    
    T& operator[](const Key &key);
    const T operator[](const Key &key) const;
    
    iterator begin();
    const_iterator begin() const;
    const_iterator constBegin() const;

    iterator find(const Key &key);
    const_iterator find(const Key &key) const;
    const_iterator constFind(const Key &key) const;

    iterator end();
    const_iterator end() const;
    const_iterator constEnd() const;
    
    int capacity() const;
    
    void clear();
    
    bool contains(const Key &key) const;
    
    int count(const Key &key) const;
    int count() const;

    bool empty() const;

    iterator insert(const Key &key, const T &value);
    iterator insertMulti(const Key &key, const T &value);

    bool isEmpty() const;
    
    const Key key(const T & value) const;
    const Key key(const T & value, const Key &defaultKey) const;
    
    QList<Key> keys() const;
    QList<Key> keys(const T &value) const;

    int remove(const Key &key);
    
    void reserve(int size);
    
    int size() const;
    
    void squeeze();
    
    T take (const Key &key);

    QList<Key> uniqueKeys() const;

    ChunkedHash<Key,T,N>& unite(const ChunkedHash<Key,T,N> &other);

    const T value(const Key &key) const;
    const T value(const Key &key, const T &defaultValue) const;
    
    QList<T> values() const;
    QList<T> values(const Key &key) const;
    
private:
    /** All of the sub-hashes in this hash */
    QVector< QHash<Key,T> > _chunks;
    
    /** The index of the chunk that contains a particular key */
    QHash<Key,qint32> key_to_chunkidx;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

////////////
//////////// Implementation of ChunkedHash::iterator
//////////// 

/** Empty constructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::iterator::iterator() 
                     : chunks(0), current_chunk(0)
{}

/** Copy constructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::iterator::iterator(
                            const typename ChunkedHash<Key,T,N>::iterator &other)
                     : chunks(other.chunks),
                       current_chunk(other.current_chunk),
                       current_it(other.current_it)
{}

/** Destructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::iterator::~iterator()
{}

/** Copy assignment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator& ChunkedHash<Key,T,N>::iterator::operator=(
                            const typename ChunkedHash<Key,T,N>::const_iterator &other)
{
    chunks = other.chunks;
    current_chunk = other.current_chunk;
    current_it = other.current_it;
    return *this;
}

/** Comparison operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::iterator::operator==(
                    const typename ChunkedHash<Key,T,N>::iterator &other) const
{
    return chunks == other.chunks and
           current_it == other.current_it;
}

/** Comparison operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::iterator::operator!=(
                    const typename ChunkedHash<Key,T,N>::iterator &other) const
{
    return not this->operator==(other);
}

/** Compare to a const_iterator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::iterator::operator==(
                    const typename ChunkedHash<Key,T,N>::const_iterator &other) const
{
    return chunks == other.chunks and
           current_it == other.current_it;
}

/** Compare to a const_iterator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::iterator::operator!=(
                    const typename ChunkedHash<Key,T,N>::const_iterator &other) const
{
    return not this->operator==(other);
}

/** Return the key of the current item */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const Key& ChunkedHash<Key,T,N>::iterator::key() const
{
    return current_it.key();
}

/** Return the value of the current item */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
T& ChunkedHash<Key,T,N>::iterator::value() const
{
    return current_it.value();
}

/** Return the value of the current item */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
T& ChunkedHash<Key,T,N>::iterator::operator*() const
{
    return *(current_it);
}

/** Return a pointer to the current value */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
T* ChunkedHash<Key,T,N>::iterator::operator->() const
{
    return current_it.operator->();
}

/** Pre-increment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator& ChunkedHash<Key,T,N>::iterator::operator++()
{
    if (chunks == 0)
        return *this;

    ++current_it;
    
    Q_ASSERT( current_chunk < chunks->count() );
    
    while (current_it == (*chunks)[current_chunk].end())
    {
        //are we on the last chunk?
        if (current_chunk >= chunks->count()-1)
        {
            //we are - we've reached the end of the hash
            return *this;
        }
        else
        {
            //move onto the next chunk
            ++current_chunk;
            current_it = (*chunks)[current_chunk].begin();
        }
    }

    return *this;
}

/** Post-increment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator ChunkedHash<Key,T,N>::iterator::operator++(int)
{
    typename ChunkedHash<Key,T,N>::iterator ret(*this);
    
    return ++ret;
}

/** Increment this iterator by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator& 
ChunkedHash<Key,T,N>::iterator::operator+=(int j)
{
    if (j < 0)
        return (*this -= -j);
        
    for (int i=0; i<j; ++i)
    {
        ++(*this);
    }
    
    return *this;
}

/** Return a copy of this iterator that has been advanced by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator 
ChunkedHash<Key,T,N>::iterator::operator+(int j) const
{
    typename ChunkedHash<Key,T,N>::iterator ret(*this);
    ret += j;
    return ret;
}

/** Pre-decrement operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator& 
ChunkedHash<Key,T,N>::iterator::operator--()
{
    if (chunks == 0)
        return *this;

    Q_ASSERT( current_chunk < chunks->count() );

    while (current_it == (*chunks)[current_chunk].begin())
    {
        //are we on the first item?
        if (current_chunk == 0)
        {
            //we are - that's the end of this hash
            return *this;
        }
        else
        {
            //move onto the previous chunk
            --current_chunk;
            
            current_it = (*chunks)[current_chunk].end();
        }
    }

    --current_it;

    return *this;
}

/** Post-decrement operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator 
ChunkedHash<Key,T,N>::iterator::operator--( int)
{
    typename ChunkedHash<Key,T,N>::iterator ret(*this);
    
    return --ret;
}

/** Decrement this iterator by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator&
ChunkedHash<Key,T,N>::iterator::operator-=( int j )
{
    if (j < 0)
        return (*this += -j);
        
    for (int i=0; i<j; ++i)
    {
        --(*this);
    }
    
    return *this;
}

/** Return a copy of this iterator that has been stepped back by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator 
ChunkedHash<Key,T,N>::iterator::operator-(int j) const
{
    typename ChunkedHash<Key,T,N>::iterator ret(*this);
    ret += j;
    return ret;
}

////////////
//////////// Implementation of ChunkedHash::const_iterator
//////////// 

/** Empty constructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::const_iterator::const_iterator() 
                     : chunks(0), current_chunk(0)
{}

/** Copy constructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::const_iterator::const_iterator(
                             const typename ChunkedHash<Key,T,N>::const_iterator &other)
                     : chunks(other.chunks),
                       current_chunk(other.current_chunk),
                       current_it(other.current_it)
{}

/** Destructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::const_iterator::~const_iterator()
{}

/** Copy assignment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator& 
ChunkedHash<Key,T,N>::const_iterator::operator=(
                            const typename ChunkedHash<Key,T,N>::iterator &other)
{
    chunks = other.chunks;
    current_chunk = other.current_chunk;
    current_it = other.current_it;
    return *this;
}

/** Copy assignment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator& 
ChunkedHash<Key,T,N>::const_iterator::operator=(
                            const typename ChunkedHash<Key,T,N>::const_iterator &other)
{
    chunks = other.chunks;
    current_chunk = other.current_chunk;
    current_it = other.current_it;
    return *this;
}

/** Comparison operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::const_iterator::operator==(
                    const typename ChunkedHash<Key,T,N>::const_iterator &other) const
{
    return chunks == other.chunks and
           current_it == other.current_it;
}

/** Comparison operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::const_iterator::operator!=(
                    const typename ChunkedHash<Key,T,N>::const_iterator &other) const
{
    return not this->operator==(other);
}

/** Compare to a iterator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::const_iterator::operator==(
                    const typename ChunkedHash<Key,T,N>::iterator &other) const
{
    return chunks == other.chunks and 
           current_it == other.current_it;
}

/** Compare to a iterator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::const_iterator::operator!=(
                    const typename ChunkedHash<Key,T,N>::iterator &other) const
{
    return not this->operator==(other);
}

/** Return the key of the current item */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const Key& ChunkedHash<Key,T,N>::const_iterator::key() const
{
    return current_it.key();
}

/** Return the value of the current item */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const T& ChunkedHash<Key,T,N>::const_iterator::value() const
{
    return current_it.value();
}

/** Return the value of the current item */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const T& ChunkedHash<Key,T,N>::const_iterator::operator*() const
{
    return *(current_it);
}

/** Return a pointer to the current value */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const T* ChunkedHash<Key,T,N>::const_iterator::operator->() const
{
    return current_it.operator->();
}

/** Pre-increment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator& 
ChunkedHash<Key,T,N>::const_iterator::operator++()
{
    if (chunks == 0)
        return *this;

    ++current_it;
    
    Q_ASSERT( current_chunk < chunks->count() );
    
    while (current_it == (*chunks)[current_chunk].end())
    {
        //are we on the last chunk?
        if (current_chunk >= chunks->count()-1)
        {
            //we are - we've reached the end of the hash
            return *this;
        }
        else
        {
            //move onto the next chunk
            ++current_chunk;
            current_it = (*chunks)[current_chunk].constBegin();
        }
    }

    return *this;
}

/** Post-increment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator 
ChunkedHash<Key,T,N>::const_iterator::operator++(int)
{
    typename ChunkedHash<Key,T,N>::const_iterator ret(*this);
    
    return ++ret;
}

/** Increment this const_iterator by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator& 
ChunkedHash<Key,T,N>::const_iterator::operator+=(int j)
{
    if (j < 0)
        return (*this -= -j);
        
    for (int i=0; i<j; ++i)
    {
        ++(*this);
    }
    
    return *this;
}

/** Return a copy of this const_iterator that has been advanced by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator 
ChunkedHash<Key,T,N>::const_iterator::operator+(int j) const
{
    typename ChunkedHash<Key,T,N>::const_iterator ret(*this);
    ret += j;
    return ret;
}

/** Pre-decrement operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator& 
ChunkedHash<Key,T,N>::const_iterator::operator--()
{
    if (chunks == 0)
        return *this;

    Q_ASSERT( current_chunk < chunks->count() );

    while (current_it == (*chunks)[current_chunk].begin())
    {
        //are we on the first item?
        if (current_chunk == 0)
        {
            //we are - that's the end of this hash
            return *this;
        }
        else
        {
            //move onto the previous chunk
            --current_chunk;
            
            current_it = (*chunks)[current_chunk].constEnd();
        }
    }

    --current_it;

    return *this;
}

/** Post-decrement operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator 
ChunkedHash<Key,T,N>::const_iterator::operator--( int)
{
    typename ChunkedHash<Key,T,N>::const_iterator ret(*this);
    
    return --ret;
}

/** Decrement this const_iterator by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator&
ChunkedHash<Key,T,N>::const_iterator::operator-=( int j )
{
    if (j < 0)
        return (*this += -j);
        
    for (int i=0; i<j; ++i)
    {
        --(*this);
    }
    
    return *this;
}

/** Return a copy of this const_iterator that has been stepped back by 'j' steps */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator 
ChunkedHash<Key,T,N>::const_iterator::operator-(int j) const
{
    typename ChunkedHash<Key,T,N>::const_iterator ret(*this);
    ret += j;
    return ret;
}

////////////
//////////// Implementation of ChunkedHash
//////////// 

/** Construct an empty hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::ChunkedHash()
{}

/** Construct from a QHash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::ChunkedHash(const QHash<Key,T> &hash)
{
    this->reserve( hash.count() );

    for (typename QHash<Key,T>::const_iterator it = hash.constBegin();
         it != hash.constEnd();
         ++it)
    {
        this->insert( it.key(), it.value() );
    }
}

/** Copy constructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::ChunkedHash(const ChunkedHash<Key,T,N> &other)
                     : _chunks(other._chunks),
                       key_to_chunkidx(other.key_to_chunkidx)
{}

/** Destructor */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>::~ChunkedHash()
{}

/** Copy assignment operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>& ChunkedHash<Key,T,N>::operator=(const ChunkedHash<Key,T,N> &other)
{
    _chunks = other._chunks;
    key_to_chunkidx = other.key_to_chunkidx;
    return *this;
}

/** Copy from a QHash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>& ChunkedHash<Key,T,N>::operator=(const QHash<Key,T> &hash)
{
    return this->operator=( ChunkedHash<Key,T,N>(hash) );
}

/** Comparison operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::operator==(const ChunkedHash<Key,T,N> &other) const
{
    return _chunks == other._chunks;
}

/** Comparison operator */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::operator!=(const ChunkedHash<Key,T,N> &other) const
{
    return not this->operator==(other);
}

/** Return a modifiable reference to the object associated with the key 'key'.
    If the hash doesn't contain any item with this key, then a default-constructed
    item is inserted into the hash and then returned */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
T& ChunkedHash<Key,T,N>::operator[](const Key &key)
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
    {
        this->insert( key, T() );
        idx = key_to_chunkidx.value(key, -1);
        
        Q_ASSERT( idx != -1 );
    }
    
    return _chunks.data()[idx][key];
}

/** Return the object associated with the key 'key', or a default constructed
    value if there is no such object in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const T ChunkedHash<Key,T,N>::operator[](const Key &key) const
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
        return T();
    else
        return _chunks.constData()[idx][key];
}

/** Return an iterator pointing to the first item in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator ChunkedHash<Key,T,N>::begin()
{
    typename ChunkedHash<Key,T,N>::iterator it;
   
    if (key_to_chunkidx.isEmpty())
        return it;
   
    it.chunks = &_chunks;
    it.current_chunk = 0;

    qint32 nchunks = _chunks.count();
    QHash<Key,T> *chunks_array = _chunks.data();

    for (qint32 i=0; i<nchunks; ++i)
    {
        if (chunks_array[i].count() != 0)
        {
            it.current_chunk = i;
            it.current_it = chunks_array[i].begin();
            break;
        }
    }
    
    return it;
}

/** Return a const_iterator pointing to the first item in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator ChunkedHash<Key,T,N>::begin() const
{
    typename ChunkedHash<Key,T,N>::const_iterator it;
   
    if (key_to_chunkidx.isEmpty())
        return it;

    it.chunks = &_chunks;
    it.current_chunk = 0;

    qint32 nchunks = _chunks.count();
    const QHash<Key,T> *chunks_array = _chunks.constData();
   
    for (qint32 i=0; i<nchunks; ++i)
    {
        if (chunks_array[i].count() != 0)
        {
            it.current_chunk = i;
            it.current_it = chunks_array[i].constBegin();
            break;
        }
    }
    
    return it;
}

/** Return a const_iterator pointing to the first item in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator ChunkedHash<Key,T,N>::constBegin() const
{
    return this->begin();
}

/** Return an iterator pointing to the most recently inserted value with key 'key',
    or end() if no such item exists */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator ChunkedHash<Key,T,N>::find(const Key &key)
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
        return this->end();
        
    typename ChunkedHash<Key,T,N>::iterator it;
    
    it.chunks = &_chunks;
    it.current_chunk = idx;
    it.current_it = _chunks[idx].find(key);
    
    return it;
}

/** Return an iterator pointing to the most recently inserted value with key 'key',
    or end() if no such item exists */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator 
ChunkedHash<Key,T,N>::find(const Key &key) const
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
        return this->end();
        
    typename ChunkedHash<Key,T,N>::const_iterator it;
    
    it.chunks = &_chunks;
    it.current_chunk = idx;
    it.current_it = _chunks[idx].find(key);
    
    return it;
}

/** Return an iterator pointing to the most recently inserted value with key 'key',
    or constEnd() if no such item exists */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator 
ChunkedHash<Key,T,N>::constFind(const Key &key) const
{
    return this->find(key);
}

/** Return an iterator that points one beyond the last item in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator ChunkedHash<Key,T,N>::end()
{
    typename ChunkedHash<Key,T,N>::iterator it;

    if (key_to_chunkidx.isEmpty())
        return it;
        
    it.chunks = &_chunks;
    it.current_chunk = _chunks.count() - 1;
    
    it.current_it = _chunks[ _chunks.count() - 1 ].end();
    
    return it;
}

/** Return an iterator that points one beyond the last item in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator ChunkedHash<Key,T,N>::end() const
{
    typename ChunkedHash<Key,T,N>::const_iterator it;

    if (key_to_chunkidx.isEmpty())
        return it;
        
    it.chunks = &_chunks;
    it.current_chunk = _chunks.count() - 1;
    
    it.current_it = _chunks[ _chunks.count() - 1 ].end();
    
    return it;
}

/** Return an iterator that points one beyond the last item in the hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::const_iterator ChunkedHash<Key,T,N>::constEnd() const
{
    return this->end();
}

/** Return the capacity of the hash (number of items that can be saved
    without requiring a new memory allocation) */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
int ChunkedHash<Key,T,N>::capacity() const
{
    int c = 0;
    
    int nchunks = _chunks.count();
    
    for (int i=0; i<nchunks; ++i)
    {
        c += _chunks.constData()[i].capacity();
    }
    
    return qMin(c, key_to_chunkidx.capacity());
}

/** Completely clear this hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
void ChunkedHash<Key,T,N>::clear()
{
    _chunks.clear();
    key_to_chunkidx.clear();
}

/** Return whether or not this has contains an item associated with the key 'key' */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::contains(const Key &key) const
{
    return key_to_chunkidx.contains(key);
}

/** Return the number of items associated with the key 'key' */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
int ChunkedHash<Key,T,N>::count(const Key &key) const
{
    return key_to_chunkidx.count(key);
}

/** Return the number of items in this hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
int ChunkedHash<Key,T,N>::count() const
{
    return key_to_chunkidx.count();
}

/** Return whether or not this hash is empty */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::empty() const
{
    return key_to_chunkidx.empty();
}

/** Insert the item with value 'value' into the hash with key 'key'. Any
    existing value with that key is overwritten */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator 
ChunkedHash<Key,T,N>::insert(const Key &key, const T &value)
{
    int idx = key_to_chunkidx.value(key, -1);

    typename ChunkedHash<Key,T,N>::iterator it;
    
    if (idx == -1)
    {
        //this is a new key
        int nchunks = _chunks.count();
        QHash<Key,T> *chunks_array = _chunks.data();
        
        //can we fit this item into an existing chunk?
        for (int i=0; i<nchunks; ++i)
        {
            if (chunks_array[i].count() < N)
            {
                key_to_chunkidx.insert(key, i);
                
                it.chunks = &(_chunks);
                it.current_chunk = i;
                it.current_it = chunks_array[i].insertMulti(key, value);
                
                return it;
            }
        }
        
        //we need to add a new chunk
        _chunks.append( QHash<Key,T>() );
        chunks_array = _chunks.data();
        chunks_array[nchunks].reserve(N);
        
        key_to_chunkidx.insert(key, nchunks);
        
        it.chunks = &(_chunks);
        it.current_chunk = nchunks;
        it.current_it = chunks_array[nchunks].insert(key, value);
    }
    else
    {
        //this will replace an existing value
        it.chunks = &(_chunks);
        it.current_chunk = idx;
        it.current_it = _chunks[idx].find(key);
        it.current_it.value() = value;
    }

    return it;
}

/** Insert a new item into the hash with key 'key' - this will create a new
    item if there is already an item with this key */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
typename ChunkedHash<Key,T,N>::iterator 
ChunkedHash<Key,T,N>::insertMulti(const Key &key, const T &value)
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
    {
        //this item doesn't exist in the hash
        return this->insert(key, value);
    }
    else
    {
        //add this value to the same chunk that contains the other value(s)
        typename ChunkedHash<Key,T,N>::iterator it;
        
        it.chunks = &_chunks;
        it.current_chunk = idx;
        it.current_it = _chunks[idx].insertMulti(key, value);
        
        key_to_chunkidx.insert(key, idx);
        
        return it;
    }
}

/** Return whether or not this hash is empty */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
bool ChunkedHash<Key,T,N>::isEmpty() const
{
    return key_to_chunkidx.isEmpty();
}

/** Return the first key associated with the value 'value', or
    'defaultKey' if there is no such key */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const Key ChunkedHash<Key,T,N>::key(const T &value, const Key &defaultKey) const
{
    int nchunks = _chunks.count();
    const QHash<Key,T> *chunks_array = _chunks.constData();
    
    for (int i=0; i<nchunks; ++i)
    {
        for (typename QHash<Key,T>::const_iterator it = chunks_array[i].constBegin();
             it != chunks_array[i].constEnd();
             ++it)
        {
            if (it.value() == value)
                return it.key();
        }
    }
    
    return defaultKey;
}

/** Return the first key associated with the value 'value' */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const Key ChunkedHash<Key,T,N>::key(const T & value) const
{
    return this->key( value, Key() );
}

/** Return all of the keys of this hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QList<Key> ChunkedHash<Key,T,N>::keys() const
{
    return key_to_chunkidx.keys();
}

/** Return all of the keys that match the value 'value' */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QList<Key> ChunkedHash<Key,T,N>::keys(const T &value) const
{
    QList<Key> ks;
    
    int nchunks = _chunks.count();
    const QHash<Key,T> *chunks_array = _chunks.constData();
    
    for (int i=0; i<nchunks; ++i)
    {
        for (typename QHash<Key,T>::const_iterator it = chunks_array[i].constBegin();
             it != chunks_array[i].constEnd();
             ++it)
        {
            if (it.value() == value)
                ks.append(it.key());
        }
    }

    return ks;
}

/** Remove all items associated with the key 'key' - return the number of items
    removed */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
int ChunkedHash<Key,T,N>::remove(const Key &key)
{
    //all values with the same key are in the same chunk
    int idx = key_to_chunkidx.value(key, -1);
    key_to_chunkidx.remove(key);
    
    if (idx != -1)
    {
        return _chunks[idx].remove(key);
    }
    else
        return 0;
}

/** Reserve at least 'size' elements - this ensures that at least
    'size' elements can be added to the has without a memory allocation */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
void ChunkedHash<Key,T,N>::reserve(int size)
{
    key_to_chunkidx.reserve(size);
    
    int nchunks = (size / N) + 1;
    
    if (_chunks.count() < nchunks)
    {
        _chunks.resize( nchunks );
        
        for (int i=0; i<nchunks; ++i)
        {
            _chunks[i].reserve(N);
        }
    }
}

/** Return the number of items in this hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
int ChunkedHash<Key,T,N>::size() const
{
    return key_to_chunkidx.size();
}

/** Squeeze this hash - this reduces its memory to the minimum
    required to hold all of the current values */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
void ChunkedHash<Key,T,N>::squeeze()
{
    //try to move items to fill up any gaps in previous chunks
}

/** Take the (most-recently inserted) item associated with the key 'key'
    from this hash and return it */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
T ChunkedHash<Key,T,N>::take(const Key &key)
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
        return T();

    key_to_chunkidx.take(key);
        
    return _chunks[idx].take(key);
}

/** Return a list of all of the unique keys in this hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QList<Key> ChunkedHash<Key,T,N>::uniqueKeys() const
{
    return key_to_chunkidx.uniqueKeys();
}

/** Unite this has with 'other' */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
ChunkedHash<Key,T,N>& ChunkedHash<Key,T,N>::unite(const ChunkedHash<Key,T,N> &other)
{
    for (typename ChunkedHash<Key,T,N>::const_iterator it = other.constBegin();
         it != other.constEnd();
         ++it)
    {
        this->insertMulti( it.key() , it.value() );
    }
    
    return *this;
}

/** Return the value associated with the key 'key', or a default constructed
    value if no such value exists */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const T ChunkedHash<Key,T,N>::value(const Key &key) const
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
        return T();
    else
        return _chunks.constData()[idx].value(key);
}

/** Return the value associated with the key 'key', or 'defaultValue'
    if no such value exists */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const T ChunkedHash<Key,T,N>::value(const Key &key, const T &defaultValue) const
{
    int idx = key_to_chunkidx.value(-1);
    
    if (idx == -1)
        return defaultValue;
    else
        return _chunks.constData()[idx].value(key);
}

/** Return the list of all values in this hash */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QList<T> ChunkedHash<Key,T,N>::values() const
{
    QList<T> vals;
    
    int nchunks = _chunks.count();
    const QHash<Key,T> *chunks_array = _chunks.constData();
    
    for (int i=0; i<nchunks; ++i)
    {
        vals += chunks_array[i].values();
    }
    
    return vals;
}

/** Return all values that are associated with the key 'key' */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QList<T> ChunkedHash<Key,T,N>::values(const Key &key) const
{
    int idx = key_to_chunkidx.value(key, -1);
    
    if (idx == -1)
        return QList<T>();
    else
        return _chunks[idx].values(key);
}

namespace detail
{

template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
const void* get_shared_container_pointer(const SireBase::ChunkedHash<Key,T,N> &hash)
{
    if (hash.empty())
        return 0;
    else
        return hash._chunks.constData();
}

template<class Key, class T, int N>
struct GetChunkedHashPointer
{
    static bool isEmpty(const ChunkedHash<Key,T,N> &hash)
    {
        return hash.empty();
    }

    static const void* value(const ChunkedHash<Key,T,N> &hash)
    {
        return get_shared_container_pointer<Key,T,N>(hash);
    }

    static void load(QDataStream &ds, ChunkedHash<Key,T,N> &hash)
    {
        ds >> hash;
    }
    
    static void save(QDataStream &ds, const ChunkedHash<Key,T,N> &hash)
    {
        ds << hash;
    }
};

} // end of namespace detail

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Serialise to a binary datastream */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireBase::ChunkedHash<Key,T,N> &hash)
{
    //this streams out the data using the same format as QHash
    ds << qint32( hash.count() );
    
    for (typename SireBase::ChunkedHash<Key,T,N>::const_iterator it = hash.constBegin();
         it != hash.constEnd();
         ++it)
    {
        ds << it.key() << it.value();
    }
    
    return ds;
}

/** Extract from a binary datastream */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireBase::ChunkedHash<Key,T,N> &hash)
{
    //this reads in using the same format as QHash
    qint32 count;
    
    ds >> count;
    
    hash.clear();
    hash.reserve(count);
    
    for (qint32 i=0; i<count; ++i)
    {
        Key key;
        T value;
        ds >> key >> value;
        
        hash.insert(key, value);
    }

    return ds;
}

/** Serialise to a binary datastream */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
SireStream::SharedDataStream& 
operator<<(SireStream::SharedDataStream &sds, const SireBase::ChunkedHash<Key,T,N> &hash)
{
    sds.sharedSaveContainer< SireBase::ChunkedHash<Key,T,N>, 
                             SireBase::detail::GetChunkedHashPointer<Key,T,N> >(hash);
                            
    return sds;
}

/** Extract from a binary datastream */
template<class Key, class T, int N>
SIRE_OUTOFLINE_TEMPLATE
SireStream::SharedDataStream& 
operator>>(SireStream::SharedDataStream &sds, SireBase::ChunkedHash<Key,T,N> &hash)
{
    sds.sharedLoadContainer< SireBase::ChunkedHash<Key,T,N>, 
                             SireBase::detail::GetChunkedHashPointer<Key,T,N> >(hash);
                            
    return sds;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
