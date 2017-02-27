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

#ifndef SIREBASE_REFCOUNTDATA_H
#define SIREBASE_REFCOUNTDATA_H

#include <tbb/atomic.h>
#include <tbb/spin_mutex.h>

#include "sireglobal.h"

#include <QDebug>

SIRE_BEGIN_HEADER

// #define SIRE_USE_REFCOUNT_MUTEX 1

namespace SireBase
{

/** This class provides a reference count that be used by objects
    that are to be held by a shared pointer, e.g. 
    SharedDataPointer, SharedPolyPointer etc.
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT RefCountData
{
public:
    RefCountData();
    ~RefCountData();
    
    class Counter
    {
    public:
        Counter();
        ~Counter();
    
        int load() const;
        int refCount() const;
        
        bool ref();
        bool deref();
        
        bool hasSingleReference() const;
        bool hasMultipleReferences() const;
        
        bool isNotReferenced() const;
        
    private:
        void doubleDereferenced() const;

        /** The atomic holding the reference count */
        tbb::atomic<int> refcount;

        #ifdef SIRE_USE_REFCOUNT_MUTEX
            /** A mutex used to protect updates to the reference count - 1 byte */
            tbb::spin_mutex mutex;
        #endif
    };

    bool operator==(const RefCountData &other) const;
    bool operator!=(const RefCountData &other) const;

    /** Actual counter - separated to follow the API of QSharedData,
        and also so that the API for RefCountData doesn't pollute any
        inheriting class */
    Counter ref;
    
private:
    RefCountData(const RefCountData &other);
    RefCountData& operator=(RefCountData &other);
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the current value of the refcount - note this may change! */
inline int RefCountData::Counter::load() const
{
    #ifdef SIRE_USE_REFCOUNT_MUTEX
        const_cast<Counter*>(this)->mutex.lock();
    #endif
    
    int value = refcount;
    
    #ifdef SIRE_USE_REFCOUNT_MUTEX
        const_cast<Counter*>(this)->mutex.unlock();
    #endif

    return value;
}

/** Return the current reference count for the object */
inline int RefCountData::Counter::refCount() const
{
    return load();
}

/** Return whether or not there are multiple references to the object */
inline bool RefCountData::Counter::hasMultipleReferences() const
{
    return load() > 1;
}

/** Return whether or not the object has only a single reference */
inline bool RefCountData::Counter::hasSingleReference() const
{
    return load() == 1;
}

/** Return whether or not the object is current unreferenced
    (has a reference count of zero) */
inline bool RefCountData::Counter::isNotReferenced() const
{
    return load() <= 0;
}

/** Increase the reference count by one */
inline bool RefCountData::Counter::ref()
{
    #ifdef SIRE_USE_REFCOUNT_MUTEX
        mutex.lock();
    #endif
    
    int oldval = refcount.fetch_and_increment();
    
    #ifdef SIRE_USE_REFCOUNT_MUTEX
        mutex.unlock();
    #endif

    return (oldval >= 0);
}

/** Decrease the reference count by one - this will raise
    an error if the reference count drops below 0 */
inline bool RefCountData::Counter::deref()
{
    #ifdef SIRE_USE_REFCOUNT_MUTEX
        mutex.lock();
    #endif
    
    int oldval = refcount.fetch_and_decrement();
    
    #ifdef SIRE_USE_REFCOUNT_MUTEX
        mutex.unlock();
    #endif
    
    if (oldval <= 0)
    {
        this->doubleDereferenced();
    }

    return (oldval > 1);
}

namespace detail
{
    tbb::spin_mutex* get_shared_null_mutex();
}

/** This function creates a single shared null-constructed instance
    of T(), which can be used as the global null for shared pointers
    of type T */
template<class T>
SIRE_INLINE_TEMPLATE
T* create_shared_null()
{
    static T *shared_null = 0;
    
    if (not shared_null)
    {
        //speculatively create the new shared_null. This is not behind
        //the mutex, as the constructed object may creata another object
        //that calls create_not_refcounted_shared_null
        T *my_null = new T();

        auto mutex = detail::get_shared_null_mutex();
        tbb::spin_mutex::scoped_lock lock(*mutex);
        
        if (not shared_null)
        {
            shared_null = my_null;
            shared_null->ref.ref();
        }
        else
        {
            lock.release();
            delete my_null;
        }
    }
    
    return shared_null;
}

template<class T>
SIRE_INLINE_TEMPLATE
T* create_not_refcounted_shared_null()
{
    static T *shared_null = 0;
    
    if (not shared_null)
    {
        //speculatively create the new shared_null. This is not behind
        //the mutex, as the constructed object may creata another object
        //that calls create_not_refcounted_shared_null        
        T *my_null = new T();

        auto mutex = detail::get_shared_null_mutex();
        tbb::spin_mutex::scoped_lock lock(*mutex);
        
        if (not shared_null)
        {
            shared_null = my_null;
        }
        else
        {
            lock.release();
            delete my_null;
        }
    }
    
    return shared_null;
}


#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
