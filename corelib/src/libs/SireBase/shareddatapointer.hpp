/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREBASE_SHAREDDATAPOINTER_HPP
#define SIREBASE_SHAREDDATAPOINTER_HPP

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This is a version of QSharedDataPointer that has been modified to
    prevent copies when being passed references to objects that
    are already held by other SharedDataPointers

    This class is copied from QSharedDataPointer in qshareddata.h
    from the Trolltech Qt 4.2.1 distribution (C) Trolltech.com

    Modifications from the original QSharedDataPointer are highlighted,
    and are (C) Christopher Woods. (obviously all QSharedDataPointer are
    changed to SharedDataPointer)

*/
template <class T>
class SIREBASE_EXPORT SharedDataPointer
{
public:

    typedef T element_type;
    typedef T value_type;
    typedef T* pointer;

    SharedDataPointer();
    ~SharedDataPointer();

    explicit SharedDataPointer(T *data);
    SharedDataPointer(const T &obj);

    SharedDataPointer(const SharedDataPointer<T> &o);
    SharedDataPointer(SharedDataPointer &&o);

    SharedDataPointer<T>& operator=(T *o);
    SharedDataPointer<T>& operator=(const T &obj);

    SharedDataPointer<T>& operator=(const SharedDataPointer<T> &o);
    SharedDataPointer<T> &operator=(SharedDataPointer<T> &&other);

    SharedDataPointer<T>& operator=(int);
    
    void detach();
    
    T& operator*();
    const T& operator*() const;
    T* operator->();
    const T* operator->() const;
    
    operator T*();
    operator const T*() const;
    
    T *data();
    const T* data() const;
    const T* constData() const;

    const T& read() const;
    T& write();

    bool operator!() const;

    bool operator==(const SharedDataPointer<T> &other) const;
    bool operator!=(const SharedDataPointer<T> &other) const;

    bool operator==(const T *other_ptr) const;
    bool operator!=(const T *other_ptr) const;

private:
    void detach_helper();

    T *d;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::SharedDataPointer() : d(0)
{}

/** Construct from a pointer to data - this will take over ownership
    of the data, and may delete it at any time. */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::SharedDataPointer(T *adata) : d(adata)
{
    if (d) 
        d->ref.ref();
}

/** Construct from a reference to an object - if this object
    is already pointed to by a SharedDataPointer then this
    will take another reference to it, otherwise this will
    copy the object (as we will assume that it is not safe
    to delete this object!) */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::SharedDataPointer(const T &obj)
{
    T *obj_ptr = const_cast<T*>(&obj);
    
    if ( obj_ptr->ref.isNotReferenced() )
    {
        //the reference count was zero - this implies that
        //this object is not held by another SharedDataPointer,
        //(it is probably on the stack) so it is not
        //safe to use this object directly - point to a clone
        //of this object.
        d = new T(obj);
        d->ref.ref();
    }
    else
    {
        //this is held by another SharedDataPointer
        obj_ptr->ref.ref();
        d = obj_ptr;
    }
}

/** Copy constructor */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::SharedDataPointer(const SharedDataPointer<T> &other)
                     : d(other.d)
{
    if (d)
        d->ref.ref();
}

/** Move constructor */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::SharedDataPointer(SharedDataPointer<T> &&other)
                     : d(other.d)
{
    other.d = 0;
}

/** Destructor */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::~SharedDataPointer()
{ 
    if (d && !d->ref.deref()) delete d;
    d = 0;
}

/** Null assignment operator - allows you to write ptr = 0 to 
    reset the pointer to null */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>& SharedDataPointer<T>::operator=(int)
{
    if (d and not d->ref.deref())
        delete d;
        
    d = 0;
    
    return *this;
}

/** Copy assigment operator - this takes over the ownership
    of 'ptr' and can delete the object at any time. */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>& SharedDataPointer<T>::operator=(T *ptr)
{
    if (ptr != d)
    {
        if (ptr)
            ptr->ref.ref();
        
        T *old = d;
        d = ptr;
        
        if (old && !old->ref.deref())
            delete old;
    }
    return *this;
}

/** Copy assignment from a reference to the object 'obj' - this
    will take a reference to this object if it is already held
    by another SharedDataPointer, or it will take a reference to
    a clone of this object */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>& SharedDataPointer<T>::operator=(const T &obj)
{
    if (d != &obj)
    {
        T *obj_ptr = const_cast<T*>(&obj);
    
        if ( obj_ptr->ref.isNotReferenced() )
        {
            //the reference count was zero - this implies that
            //this object is not held by another SharedDataPointer,
            //(it is probably on the stack) so it is not
            //safe to use this object directly - point to a clone
            //of this object.
            obj_ptr = new T(obj);

            obj_ptr->ref.ref();
            
            T *old = d;
            d = obj_ptr;
            
            if (old && !old->ref.deref())
                delete old;
        }
        else
        {
            //this is held by another SharedDataPointer
            if (&obj != d)
            {
                const_cast<T*>(&obj)->ref.ref();
                
                T *old = d;
                d = const_cast<T*>(&obj);
                
                if (old && !old->ref.deref())
                    delete old;
            }
        }
    }
    
    return *this;
}

/** Copy assignment from another SharedDataPointer<T> */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>& SharedDataPointer<T>::operator=(const SharedDataPointer<T> &other)
{
    if (other.d != d)
    {
        if (other.d)
            other.d->ref.ref();
        
        T *old = d;
        d = other.d;
        
        if (old && !old->ref.deref())
            delete old;
    }
    
    return *this;
}

/** Move assignment operator */
template<class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>& SharedDataPointer<T>::operator=(SharedDataPointer<T> &&other)
{
    qSwap(d, other.d);
    return *this;
}

/** Helper for the detach function */
template<class T>
Q_OUTOFLINE_TEMPLATE
void SharedDataPointer<T>::detach_helper()
{
    T *x = new T(*d);
    x->ref.ref();
    if (!d->ref.deref())
        delete d;
    d = x;
}

/** Detach the object pointed to by this pointer from shared storage */
template <class T>
Q_INLINE_TEMPLATE
void SharedDataPointer<T>::detach() 
{
    if (d && d->ref.hasMultipleReferences()) detach_helper();
}

/** Dereference this pointer */
template <class T>
Q_INLINE_TEMPLATE
T& SharedDataPointer<T>::operator*() 
{ 
    detach(); 
    return *d; 
}

/** Dereference this pointer */
template <class T>
Q_INLINE_TEMPLATE
const T& SharedDataPointer<T>::operator*() const 
{
    return *d; 
}

/** Dereference this pointer for reading (const-access) */
template<class T>
Q_INLINE_TEMPLATE
const T& SharedDataPointer<T>::read() const
{
    return *d;
}

/** Dereference this pointer for writing (non-const-access) */
template<class T>
Q_INLINE_TEMPLATE
T& SharedDataPointer<T>::write()
{
    detach();
    return *d;
}

/** Pointer dereference */
template <class T>
Q_INLINE_TEMPLATE
T* SharedDataPointer<T>::operator->() 
{
    detach(); 
    return d; 
}

/** Pointer dereference */
template <class T>
Q_INLINE_TEMPLATE
const T* SharedDataPointer<T>::operator->() const 
{
    return d;
}

/** Cast back to a normal pointer */
template <class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::operator T*() 
{
    detach(); 
    return d; 
}

/** Cast back to a normal pointer */
template <class T>
Q_INLINE_TEMPLATE
SharedDataPointer<T>::operator const T*() const 
{
    return d;
}

/** Return the pointer held by this shared pointer */
template <class T>
Q_INLINE_TEMPLATE
T* SharedDataPointer<T>::data()
{
    detach(); 
    return d; 
}

/** Return the pointer held by this shared pointer */
template <class T>
Q_INLINE_TEMPLATE
const T* SharedDataPointer<T>::data() const
{
    return d;
}

/** Return the pointer held by this shared pointer */
template <class T>
Q_INLINE_TEMPLATE
const T* SharedDataPointer<T>::constData() const
{
    return d;
}

/** Used to implement if (!ptr){ ... } */
template <class T>
Q_INLINE_TEMPLATE
bool SharedDataPointer<T>::operator!() const 
{
    return !d;
}

/** Comparison operator */
template <class T>
Q_INLINE_TEMPLATE
bool SharedDataPointer<T>::operator==(const SharedDataPointer<T> &other) const
{
    return d == other.d;
}

/** Comparison operator */
template <class T>
Q_INLINE_TEMPLATE
bool SharedDataPointer<T>::operator!=(const SharedDataPointer<T> &other) const
{
    return d != other.d;
}

/** Comparison operator */
template <class T>
Q_INLINE_TEMPLATE
bool SharedDataPointer<T>::operator==(const T *other_ptr) const
{
    return d == other_ptr;
}

/** Comparison operator */
template <class T>
Q_INLINE_TEMPLATE
bool SharedDataPointer<T>::operator!=(const T *other_ptr) const
{
    return d != other_ptr;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireBase

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
inline const T* get_pointer(SireBase::SharedDataPointer<T> const & p)
{
    return p.constData();
}

template<class T>
inline T* get_pointer(SireBase::SharedDataPointer<T> &p)
{
    return p.data();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
