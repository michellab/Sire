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

SIRE_BEGIN_HEADER

#define SIRE_USE_REFCOUNT_MUTEX 1

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
        auto mutex = detail::get_shared_null_mutex();
        
        tbb::spin_mutex::scoped_lock lock(*mutex);
        
        T *my_null = new T();
        
        if (not shared_null)
        {
            shared_null = my_null;
            shared_null->ref.ref();
        }
        else
            delete my_null;
    }
    
    return shared_null;
}

}

SIRE_END_HEADER

#endif
