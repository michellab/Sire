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

#include <atomic>
#include <tbb/spin_mutex.h>

#include "sireglobal.h"

#include <QDebug>

SIRE_BEGIN_HEADER

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

        void reset();

        bool hasSingleReference() const;
        bool hasMultipleReferences() const;

        bool isNotReferenced() const;

    private:
        /** The atomic holding the reference count */
        std::atomic<int> refcount;
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

/** Reset the counter so that it is set back equal to zero */
SIRE_ALWAYS_INLINE void RefCountData::Counter::reset()
{
    refcount.exchange(0);
}

/** Return the current value of the refcount - note this may change! */
SIRE_ALWAYS_INLINE int RefCountData::Counter::load() const
{
    int value = refcount;

    return value;
}

/** Return the current reference count for the object */
SIRE_ALWAYS_INLINE int RefCountData::Counter::refCount() const
{
    return load();
}

/** Return whether or not there are multiple references to the object */
SIRE_ALWAYS_INLINE bool RefCountData::Counter::hasMultipleReferences() const
{
    return load() > 1;
}

/** Return whether or not the object has only a single reference */
SIRE_ALWAYS_INLINE bool RefCountData::Counter::hasSingleReference() const
{
    return load() == 1;
}

/** Return whether or not the object is current unreferenced
    (has a reference count of zero) */
SIRE_ALWAYS_INLINE bool RefCountData::Counter::isNotReferenced() const
{
    return load() <= 0;
}

/** Increase the reference count by one */
SIRE_ALWAYS_INLINE bool RefCountData::Counter::ref()
{
    int oldval = refcount.fetch_add(1);

    return (oldval >= 0);
}

/** Decrease the reference count by one - this will raise
    an error if the reference count drops below 0 */
SIRE_ALWAYS_INLINE bool RefCountData::Counter::deref()
{
    int oldval = refcount.fetch_sub(1);

    return (oldval > 1);
}

namespace detail
{
    SIREBASE_EXPORT tbb::spin_mutex* get_shared_null_mutex();
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
