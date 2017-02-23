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

#include "refcountdata.h"

#include <QDebug>

using namespace SireBase;

/** Constructor */
RefCountData::RefCountData()
{}

/** Copying is forbidden */
RefCountData::RefCountData(const RefCountData&o)
{}

/** Destructor */
RefCountData::~RefCountData()
{}

/** Copy assignment is forbidden */
RefCountData& RefCountData::operator=(RefCountData&)
{
    return *this;
}

/** Only the same if exactly the same object */
bool RefCountData::operator==(const RefCountData &other) const
{
    return this == &other;
}

/** Only the same if exactly the same object */
bool RefCountData::operator!=(const RefCountData &other) const
{
    return not operator==(other);
}

/** Return the current value of the refcount - note this may change! */
int RefCountData::Counter::load() const
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

/** Constructor */
RefCountData::Counter::Counter()
{
    refcount = 0;
}

/** Destructor */
RefCountData::Counter::~Counter()
{}

/** Return the current reference count for the object */
int RefCountData::Counter::refCount() const
{
    return load();
}

/** Return whether or not there are multiple references to the object */
bool RefCountData::Counter::hasMultipleReferences() const
{
    return load() > 1;
}

/** Return whether or not the object has only a single reference */
bool RefCountData::Counter::hasSingleReference() const
{
    return load() == 1;
}

/** Return whether or not the object is current unreferenced
    (has a reference count of zero) */
bool RefCountData::Counter::isNotReferenced() const
{
    return load() <= 0;
}

/** Increase the reference count by one */
bool RefCountData::Counter::ref()
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
bool RefCountData::Counter::deref()
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
        //we just reduced the count to zero a second time!
        qDebug() << "WARNING - PROGRAM BUG - REDUCED REFCOUNT TO 0 TWICE!";
    }
    
    return (oldval > 0);
}
