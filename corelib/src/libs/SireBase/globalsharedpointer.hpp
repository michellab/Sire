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

#ifndef SIREBASE_GLOBALSHAREDPOINTER_HPP
#define SIREBASE_GLOBALSHAREDPOINTER_HPP

#include <QMutex>
#include <QMutexLocker>
#include <QSet>

#include "sharedpolypointer.hpp"

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class GlobalSharedPointer;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::GlobalSharedPointer<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::GlobalSharedPointer<T>&);

namespace SireBase
{

/** This is the base class of global shared pointers */
class SIREBASE_EXPORT GlobalSharedPointerBase
{
public:
    ~GlobalSharedPointerBase();

protected:
    GlobalSharedPointerBase();

    static QMutex* registryMutex();
    static QSet<const void*>& getRegistry(const char *typname); 

    template<class T>
    static T* registerObject(const T *obj_ptr);
    
    template<class T>
    static void unregisterObject(const T *obj_ptr);
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* GlobalSharedPointerBase::registerObject(const T *obj_ptr)
{
    if (not obj_ptr)
        return const_cast<T*>(obj_ptr);

    QMutexLocker lkr( GlobalSharedPointerBase::registryMutex() );
    
    QSet<const void*> &registry = GlobalSharedPointerBase::getRegistry(
                                        SharedPolyPointerHelper<T>::typeName() );
                                        
    for (typename QSet<const void*>::const_iterator it = registry.constBegin();
         it != registry.constEnd();
         ++it)
    {
        const T *global_obj = (const T*)(*it);
        
        if ( obj_ptr == global_obj)
            //this is the same pointer
            return const_cast<T*>(obj_ptr);

        if ( SharedPolyPointerHelper<T>::equal(*obj_ptr, *global_obj) )
        {
            //the objects are equal - return the global copy
            return const_cast<T*>(global_obj);
        }
    }
    
    //there is no global object with this value
    registry.insert( obj_ptr );
    
    return const_cast<T*>(obj_ptr);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void GlobalSharedPointerBase::unregisterObject(const T *obj_ptr)
{
    if (not obj_ptr)
        return;

    QMutexLocker lkr( GlobalSharedPointerBase::registryMutex() );
    
    QSet<const void*> &registry = GlobalSharedPointerBase::getRegistry(
                                        SharedPolyPointerHelper<T>::typeName() );

    if (registry.contains(obj_ptr))
    {
        //we can only remove this from the registry if this is a unique
        //pointer
        if (obj_ptr->ref.hasSingleReference())
            registry.remove(obj_ptr);
    }
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

/** This is a SharedPolyPointer that uses a global registry
    to ensure that there is just one shared copy of an object.
    
    Note that the this can only hold a const pointer - it does
    not allow editing of the object. 
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT GlobalSharedPointer
        : private SharedPolyPointer<T>, private GlobalSharedPointerBase
{

friend QDataStream& ::operator<<<>(QDataStream&, const GlobalSharedPointer<T>&);
friend QDataStream& ::operator>><>(QDataStream&, GlobalSharedPointer<T>&);

public:

    typedef T element_type;
    typedef T value_type;
    typedef T* pointer;

    GlobalSharedPointer();
    ~GlobalSharedPointer();

    explicit GlobalSharedPointer(T *data);
    GlobalSharedPointer(const T &obj);

    GlobalSharedPointer(const SharedPolyPointer<T> &o);

    GlobalSharedPointer(const GlobalSharedPointer<T> &o);

    template<class S>
    GlobalSharedPointer(const GlobalSharedPointer<S> &o);

    template<class S>
    explicit GlobalSharedPointer(S *data);

    template<class S>
    GlobalSharedPointer(const S &obj);

    GlobalSharedPointer<T>& operator=(T *o);
    GlobalSharedPointer<T>& operator=(const T &obj);

    GlobalSharedPointer<T>& operator=(const GlobalSharedPointer<T> &o);

    template<class S>
    GlobalSharedPointer<T>& operator=(const GlobalSharedPointer<S> &o);

    GlobalSharedPointer<T>& operator=(int);

    template<class S>
    GlobalSharedPointer<T>& operator=(const S &obj);

    template<class S>
    GlobalSharedPointer<T>& operator=(S *obj);
    
    bool unique() const;
    
    const T& operator*() const;
    const T* operator->() const;
    
    operator const T*() const;
    
    const T* data() const;
    const T* constData() const;

    bool operator!() const;

    bool operator==(const GlobalSharedPointer<T> &other) const;
    bool operator!=(const GlobalSharedPointer<T> &other) const;

    bool operator==(const T *other_ptr) const;
    bool operator!=(const T *other_ptr) const;

    const char* what() const;

    template<class S>
    bool isA() const;

    template<class S>
    const S& asA() const;

private:
    void registerGlobalPointer(bool auto_detach=true);
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer()
                       : SharedPolyPointer<T>(), GlobalSharedPointerBase()
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::~GlobalSharedPointer()
{
    //are we the only pointer to this object?
    if (SharedPolyPointer<T>::unique())
    {
        //unregister this object
        GlobalSharedPointerBase::unregisterObject<T>( this->constData() );
    }
}

/** Internal function called to register the global pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void GlobalSharedPointer<T>::registerGlobalPointer(bool auto_detach)
{
    if (this->constData())
    {   
        if (auto_detach)
        {
            if (not this->unique())
                this->detach();
        }

        SharedPolyPointer<T>::operator=(
                    GlobalSharedPointerBase::registerObject<T>(this->constData()) );
    }
}

/** Construct a pointer to the object 'data' - this takes over ownership
    of the object pointed to by 'data' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(T *data)
                       : SharedPolyPointer<T>(data), GlobalSharedPointerBase()
{
    this->registerGlobalPointer();
}

/** Construct a pointer to the object 'obj' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(const T &obj)
                       : SharedPolyPointer<T>(obj), GlobalSharedPointerBase()
{
    this->registerGlobalPointer();
}

/** Construct a pointer to the object pointed to by 'o' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(const SharedPolyPointer<T> &o)
                       : SharedPolyPointer<T>(o), GlobalSharedPointerBase()
{
    this->registerGlobalPointer();
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(const GlobalSharedPointer<T> &o)
               : SharedPolyPointer<T>( static_cast<const SharedPolyPointer<T>&>(o) ), 
                 GlobalSharedPointerBase()
{
    //by definition, 'o' is already global
}

/** Casting constructor */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(const GlobalSharedPointer<S> &o)
                       : SharedPolyPointer<T>(o), GlobalSharedPointerBase()
{
    this->registerGlobalPointer();
}

/** Copy from the passed pointer - this takes over ownership of the pointer */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(S *data)
                       : SharedPolyPointer<T>(data), GlobalSharedPointerBase()
{
    this->registerGlobalPointer();
}

/** Construct from the passed object */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::GlobalSharedPointer(const S &obj)
                       : SharedPolyPointer<T>(obj), GlobalSharedPointerBase()
{
    this->registerGlobalPointer();
}

/** Copy assignment from the passed pointer - this takes over ownership of the pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(T *o)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );

    SharedPolyPointer<T>::operator=(o);
    this->registerGlobalPointer();
    
    return *this;
}

/** Copy assignment from the passed object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(const T &obj)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );

    SharedPolyPointer<T>::operator=(obj);
    this->registerGlobalPointer();
    
    return *this;
}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(
                                                const GlobalSharedPointer<T> &o)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );

    SharedPolyPointer<T>::operator=( static_cast<const SharedPolyPointer<T>&>(o) );
    
    return *this;
}

/** Copy assignment operator */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(
                                                const GlobalSharedPointer<S> &o)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );

    SharedPolyPointer<T>::operator=(o);
    this->registerGlobalPointer();

    return *this;
}

/** This is used to allow "ptr = 0" to reset the pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(int val)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );

    SharedPolyPointer<T>::operator=(val);
    
    return *this;
}

/** Copy to point at 'obj' */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(const S &obj)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );
    
    SharedPolyPointer<T>::operator=(obj);
    this->registerGlobalPointer();
    
    return *this;
}

/** Copy to point at 'obj' - this takes over ownership of the object */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>& GlobalSharedPointer<T>::operator=(S *obj)
{
    if (SharedPolyPointer<T>::unique())
        GlobalSharedPointerBase::unregisterObject( this->constData() );

    SharedPolyPointer<T>::operator=(obj);
    this->registerGlobalPointer();
    
    return *this;
}

/** Return whether or not this object is unique (and thus globally unique) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::unique() const
{
    return SharedPolyPointer<T>::unique();
}

/** Return a const reference to the object - don't dereference a null pointer! */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& GlobalSharedPointer<T>::operator*() const
{
    return SharedPolyPointer<T>::operator*();
}

/** Pointer operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* GlobalSharedPointer<T>::operator->() const
{
    return SharedPolyPointer<T>::operator->();
}

/** Cast to a raw pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalSharedPointer<T>::operator const T*() const
{
    return SharedPolyPointer<T>::constData();
}

/** Return the raw pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* GlobalSharedPointer<T>::data() const
{
    return SharedPolyPointer<T>::constData();
}

/** Return the raw pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* GlobalSharedPointer<T>::constData() const
{
    return SharedPolyPointer<T>::constData();
}

/** Used to implement if (!ptr){ ... } */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::operator!() const
{
    return SharedPolyPointer<T>::operator!();
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::operator==(const GlobalSharedPointer<T> &other) const
{
    //this is easy, as only one copy of each global object
    return this->constData() == other.constData();
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::operator!=(const GlobalSharedPointer<T> &other) const
{
    //this is easy, as only one copy of each global object
    return this->constData() != other.constData();
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::operator==(const T *other_ptr) const
{
    return SharedPolyPointer<T>::operator==(other_ptr);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::operator!=(const T *other_ptr) const
{
    return SharedPolyPointer<T>::operator!=(other_ptr);
}

/** Return the class name of the object pointed to by this pointer */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* GlobalSharedPointer<T>::what() const
{
    return SharedPolyPointer<T>::what();
}

/** Return whether or not this object is of type 'S' */
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalSharedPointer<T>::isA() const
{
    return SharedPolyPointer<T>::template isA<S>();
}

/** Return the object pointed to by this pointer cast as
    an object of type 'S'
    
    \throw SireError::invalid_cast
*/
template<class T>
template<class S>
SIRE_OUTOFLINE_TEMPLATE
const S& GlobalSharedPointer<T>::asA() const
{
    return SharedPolyPointer<T>::template asA<S>();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
QDataStream& operator<<(QDataStream &ds, const SireBase::GlobalSharedPointer<T> &obj)
{
    ds << static_cast<const SireBase::SharedPolyPointer<T>&>(obj);
    return ds;
}

template<class T>
QDataStream& operator>>(QDataStream &ds, SireBase::GlobalSharedPointer<T> &obj)
{
    obj = 0;
    
    ds >> static_cast<SireBase::SharedPolyPointer<T>&>(obj);
    
    obj.registerGlobalPointer(false);
    
    return ds;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
