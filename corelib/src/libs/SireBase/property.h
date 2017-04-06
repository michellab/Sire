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

#ifndef SIREBASE_PROPERTY_H
#define SIREBASE_PROPERTY_H

#include "sharedpolypointer.hpp"
#include "globalsharedpointer.hpp"

#include "refcountdata.h"

#include <QVariant>
#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireBase
{
class Property;
class VariantProperty;
class NullProperty;

class PropPtrBase;
class GlobalPropPtrBase;

template<class T>
class PropPtr;

template<class T>
class GlobalPropPtr;

}

QDataStream& operator<<(QDataStream&, const SireBase::Property&);
QDataStream& operator>>(QDataStream&, SireBase::Property&);

QDataStream& operator<<(QDataStream&, const SireBase::VariantProperty&);
QDataStream& operator>>(QDataStream&, SireBase::VariantProperty&);

QDataStream& operator<<(QDataStream&, const SireBase::NullProperty&);
QDataStream& operator>>(QDataStream&, SireBase::NullProperty&);

QDataStream& operator<<(QDataStream&, const SireBase::PropPtrBase&);
QDataStream& operator>>(QDataStream&, SireBase::PropPtrBase&);

QDataStream& operator<<(QDataStream&, const SireBase::GlobalPropPtrBase&);
QDataStream& operator>>(QDataStream&, SireBase::GlobalPropPtrBase&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::PropPtr<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::PropPtr<T>&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::GlobalPropPtr<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::GlobalPropPtr<T>&);

class QMutex;

namespace SireBase
{

QMutex* globalLock();

/** This is the base class of all properties that may be attached to a
    molecule. Properties are used to assign extra information to a molecule,
    which may then be carried by the molecule throughout its passage
    through the simulation. Examples of properties may include the file
    from which the molecule was read, the charge parameters on the atoms,
    the PDB code etc.

    Properties form a polymorphic hierarchy which are implicitly shared
    via SireBase::SharedPolyPointer smart pointers.

    @author Christopher Woods
*/
class SIREBASE_EXPORT Property : public RefCountData
{

friend QDataStream& ::operator<<(QDataStream&, const Property&);
friend QDataStream& ::operator>>(QDataStream&, Property&);

public:
    typedef Property ROOT;

    Property();

    Property(const Property &other);

    virtual ~Property();

    virtual Property* clone() const=0;
    virtual Property* create() const=0;
    
    virtual QString toString() const;
    
    virtual void copy(const Property &other)=0;
    virtual bool equals(const Property &other) const=0;

    virtual void save(QDataStream &ds) const=0;
    virtual void load(QDataStream &ds)=0;

    //virtual void save(XMLStream &xs) const=0;
    //virtual void load(XMLStream &xs)=0;

    virtual const char* what() const=0;

    static const char* typeName()
    {
        return "SireBase::Property";
    }

    template<class T>
    bool isA() const;

    template<class T>
    const T& asA() const;
    
    template<class T>
    T& asA();

    static const NullProperty& null();

protected:
    Property& operator=(const Property &other);

    bool operator==(const Property &other) const;
    bool operator!=(const Property &other) const;

    void throwInvalidCast(const Property &other) const;
    void throwInvalidCast(const char *typenam) const;
};

/** This is the second-to-top class of all Properties. Any instantiatable
    properties must inherit from this template class so that all of the
    virtual functions are supplied correctly, e.g. if you have a
    class called Derived, that inherits from Base, and Base derives
    from Property, then you need to inherit Derived from
    ConcreteProperty<Derived,Base>

    @author Christopher Woods
*/
template<class Derived, class Base>
class SIREBASE_EXPORT ConcreteProperty : public Base
{
public:
    ConcreteProperty();

    template<class T0>
    ConcreteProperty(const T0 &t0);

    template<class T0, class T1>
    ConcreteProperty(const T0 &t0, const T1 &t1);

    template<class T0, class T1, class T2>
    ConcreteProperty(const T0 &t0, const T1 &t1, const T2 &t2);

    template<class T0, class T1, class T2, class T3>
    ConcreteProperty(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3);

    template<class T0, class T1, class T2, class T3, class T4>
    ConcreteProperty(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3,
                     const T4 &t4);

    template<class T0, class T1, class T2, class T3, class T4, class T5>
    ConcreteProperty(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3,
                     const T4 &t4, const T5 &t5);

    template<class T0, class T1, class T2, class T3, class T4, 
             class T5, class T6>
    ConcreteProperty(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3,
                     const T4 &t4, const T5 &t5, const T6 &t6);

    template<class T0, class T1, class T2, class T3, class T4, 
             class T5, class T6, class T7>
    ConcreteProperty(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3,
                     const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7);

    virtual ~ConcreteProperty();

    ConcreteProperty<Derived,Base>& operator=(const Property &other);

    bool operator==(const Property &other) const;

    bool operator!=(const Property &other) const;

    const char* what() const;

    ConcreteProperty<Derived,Base>* clone() const;

    ConcreteProperty<Derived,Base>* create() const;

    void copy(const Property &other);

    bool equals(const Property &other) const;

    void save(QDataStream &ds) const;
    void load(QDataStream &ds);

protected:
    ConcreteProperty<Derived,Base>&
    operator=(const ConcreteProperty<Derived,Base> &other);

    bool operator==(const ConcreteProperty<Derived,Base> &other) const;
    bool operator!=(const ConcreteProperty<Derived,Base> &other) const;
};

/** This is a simple property that holds any value as a QVariant. This
    is designed to be used for metadata that doesn't need any tight
    checking (e.g. the author of the molecule file, the source of
    the coordinates, the 'header' lines etc.)

    @author Christopher Woods
*/
class SIREBASE_EXPORT VariantProperty
          : public ConcreteProperty<VariantProperty,Property>, public QVariant
{
public:
    VariantProperty();

    VariantProperty(const QVariant &value);

    VariantProperty(const Property &other);

    VariantProperty(const VariantProperty &other);

    virtual ~VariantProperty();

    VariantProperty& operator=(const QVariant &value);
    VariantProperty& operator=(const VariantProperty &other);

    bool operator==(const VariantProperty &other) const;
    bool operator!=(const VariantProperty &other) const;

    static const char* typeName();
    
    QString toString() const;
            
    template<class T>
    T convertTo() const
    {
        if (not this->canConvert<T>())
            this->throwInvalidCast( QMetaType::typeName( qMetaTypeId<T>() ) );
            
        return this->value<T>();
    }

private:
    void throwInvalidCast(const QString &typname) const;
};

/** This is a null property */
class SIREBASE_EXPORT NullProperty
              : public ConcreteProperty<NullProperty,Property>
{
public:
    NullProperty();
    
    NullProperty(const NullProperty &other);

    ~NullProperty();

    static const char* typeName();
    
    QString toString() const;
};

/** This is base class of the polymorphic pointer holder for the entire
    Property class hierarchy. This can hold implicitly 
    shared pointers to any property class.

    @author Christopher Woods
*/
class SIREBASE_EXPORT PropPtrBase
{

friend QDataStream& ::operator<<(QDataStream&, const PropPtrBase&);
friend QDataStream& ::operator>>(QDataStream&, PropPtrBase&);

public:
    ~PropPtrBase();

    bool operator==(const Property &property) const;
    bool operator!=(const Property &property) const;

    bool operator==(const PropPtrBase &other) const;
    bool operator!=(const PropPtrBase &other) const;

    void detach();

    bool unique() const;

    operator const Property&() const;

protected:
    PropPtrBase(const Property &property);
    PropPtrBase(Property *property);
    
    PropPtrBase(const PropPtrBase &other);

    PropPtrBase& operator=(const PropPtrBase &other);

    const Property& read() const;
    
    Property& write();
    Property& edit();

    static void throwCastingError(const char *got_type, const char *want_type);

private:
    /** Shared pointer to the actual property */
    SharedPolyPointer<Property> ptr;
};

/** This is base class of the global polymorphic pointer holder for the entire
    Property class hierarchy. This can hold implicitly 
    shared pointers to any property class.

    @author Christopher Woods
*/
class SIREBASE_EXPORT GlobalPropPtrBase
{

friend QDataStream& ::operator<<(QDataStream&, const GlobalPropPtrBase&);
friend QDataStream& ::operator>>(QDataStream&, GlobalPropPtrBase&);

public:
    ~GlobalPropPtrBase();

    bool operator==(const Property &property) const;
    bool operator!=(const Property &property) const;

    bool operator==(const GlobalPropPtrBase &other) const;
    bool operator!=(const GlobalPropPtrBase &other) const;

    bool unique() const;

    operator const Property&() const;

protected:
    GlobalPropPtrBase(const Property &property);
    GlobalPropPtrBase(Property *property);
    
    GlobalPropPtrBase(const GlobalPropPtrBase &other);

    GlobalPropPtrBase& operator=(const GlobalPropPtrBase &other);

    const Property& read() const;

    static void throwCastingError(const char *got_type, const char *want_type);

private:
    /** Global shared pointer to the actual property */
    GlobalSharedPointer<Property> ptr;
};

/** This is the specialised pointer that is used to hold a hierarchy
    of properties that are derived from type 'T'
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT PropPtr : public PropPtrBase
{

friend QDataStream& ::operator<<<>(QDataStream&, const PropPtr<T>&);
friend QDataStream& ::operator>><>(QDataStream&, PropPtr<T>&);

public:
    PropPtr();

    PropPtr(const T &obj);
    PropPtr(T *obj);
    
    PropPtr(const Property &property);

    PropPtr(const PropPtrBase &other);

    PropPtr(const PropPtr<T> &other);
    
    ~PropPtr();
    
    PropPtr<T>& operator=(const T &obj);
    PropPtr<T>& operator=(T *obj);
    
    PropPtr<T>& operator=(const Property &property);

    PropPtr<T>& operator=(const PropPtr<T> &other);
    PropPtr<T>& operator=(const PropPtrBase &property);

    const T* operator->() const;
    const T& operator*() const;
    
    const T& read() const;
    T& edit();

    const T* data() const;
    const T* constData() const;

    T* data();
    
    operator const T&() const;

    bool isNull() const;

    static PropPtr<T> null();

protected:
    void assertSane() const;
};

/** This is the specialised global pointer that is used to hold a hierarchy
    of properties that are derived from type 'T'
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT GlobalPropPtr : public GlobalPropPtrBase
{

friend QDataStream& ::operator<<<>(QDataStream&, const GlobalPropPtr<T>&);
friend QDataStream& ::operator>><>(QDataStream&, GlobalPropPtr<T>&);

public:
    GlobalPropPtr();

    GlobalPropPtr(const T &obj);
    GlobalPropPtr(T *obj);
    
    GlobalPropPtr(const Property &property);

    GlobalPropPtr(const GlobalPropPtrBase &other);

    GlobalPropPtr(const GlobalPropPtr<T> &other);
    
    ~GlobalPropPtr();
    
    GlobalPropPtr<T>& operator=(const T &obj);
    GlobalPropPtr<T>& operator=(T *obj);
    
    GlobalPropPtr<T>& operator=(const Property &property);

    GlobalPropPtr<T>& operator=(const GlobalPropPtr<T> &other);
    GlobalPropPtr<T>& operator=(const GlobalPropPtrBase &property);

    const T* operator->() const;
    const T& operator*() const;
    
    const T& read() const;

    const T* data() const;
    const T* constData() const;
    
    operator const T&() const;

    bool isNull() const;

    static GlobalPropPtr<T> null();

protected:
    void assertSane() const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

///////
/////// Implementation of Property template functions
///////

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Property::isA() const
{
    if (dynamic_cast<const T*>(this) != 0)
    {
        //this isA<T>() !
        return true;
    }
    else
    {
        return QLatin1String(T::typeName()) == QLatin1String(this->what());
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Property::asA() const
{
    const T* as_t = dynamic_cast<const T*>(this);
    
    if (not as_t)
    {
        if (QLatin1String(T::typeName()) == QLatin1String(this->what()))
        {
            //these are the same object - perhaps the typeinfo objects
            //are in different shared libraries?
            as_t = (const T*)(this);
        }
        else
            throwInvalidCast(T::typeName());
    }
    
    return *as_t;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& Property::asA()
{
    T* as_t = dynamic_cast<T*>(this);
    
    if (not as_t)
    {
        if (QLatin1String(T::typeName()) == QLatin1String(this->what()))
        {
            //these are the same object - perhaps the typeinfo objects
            //are in different shared libraries?
            as_t = (T*)(this);
        }
        else
            throwInvalidCast(T::typeName());
    }
       
    return *as_t;
}

///////
/////// Implementation of ConcreteProperty<T>
///////

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty() : Base()
{}

template<class Derived, class Base>
template<class T0>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0) : Base(t0)
{}

template<class Derived, class Base>
template<class T0, class T1>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1) 
                               : Base(t0,t1)
{}

template<class Derived, class Base>
template<class T0, class T1, class T2>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1,
                                                 const T2 &t2) : Base(t0,t1,t2)
{}

template<class Derived, class Base>
template<class T0, class T1, class T2, class T3>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1,
                                                 const T2 &t2, const T3 &t3) 
                               : Base(t0,t1,t2,t3)
{}

template<class Derived, class Base>
template<class T0, class T1, class T2, class T3,
         class T4>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1,
                                                 const T2 &t2, const T3 &t3,
                                                 const T4 &t4) 
                               : Base(t0,t1,t2,t3,t4)
{}

template<class Derived, class Base>
template<class T0, class T1, class T2, class T3,
         class T4, class T5>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1,
                                                 const T2 &t2, const T3 &t3,
                                                 const T4 &t4, const T5 &t5) 
                               : Base(t0,t1,t2,t3,t4,t5)
{}

template<class Derived, class Base>
template<class T0, class T1, class T2, class T3,
         class T4, class T5, class T6>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1,
                                                 const T2 &t2, const T3 &t3,
                                                 const T4 &t4, const T5 &t5,
                                                 const T6 &t6) 
                               : Base(t0,t1,t2,t3,t4,t5,t6)
{}

template<class Derived, class Base>
template<class T0, class T1, class T2, class T3,
         class T4, class T5, class T6, class T7>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::ConcreteProperty(const T0 &t0, const T1 &t1,
                                                 const T2 &t2, const T3 &t3,
                                                 const T4 &t4, const T5 &t5,
                                                 const T6 &t6, const T7 &t7) 
                               : Base(t0,t1,t2,t3,t4,t5,t6,t7)
{}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>::~ConcreteProperty()
{}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>& 
ConcreteProperty<Derived,Base>::operator=(const Property &other)
{
    const Derived* other_t = dynamic_cast<const Derived*>(&other);

    if (!other_t)
        Base::throwInvalidCast(other);

    return static_cast<Derived*>(this)->operator=(*other_t);
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
bool ConcreteProperty<Derived,Base>::operator==(const Property &other) const
{
    const Derived* other_t = dynamic_cast<const Derived*>(&other);

    if (other_t)
        return static_cast<const Derived*>(this)->operator==(*other_t);
    else
        return false;
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
bool ConcreteProperty<Derived,Base>::operator!=(const Property &other) const
{
    const Derived* other_t = dynamic_cast<const Derived*>(&other);

    if (other_t)
        return static_cast<const Derived*>(this)->operator!=(*other_t);
    else
        return true;
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
const char* ConcreteProperty<Derived,Base>::what() const
{
    return Derived::typeName();
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>* ConcreteProperty<Derived,Base>::clone() const
{
    return new Derived( static_cast<const Derived&>(*this) );
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>* ConcreteProperty<Derived,Base>::create() const
{
    return new Derived();
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
void ConcreteProperty<Derived,Base>::copy(const Property &other)
{
    ConcreteProperty<Derived,Base>::operator=(other);
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
bool ConcreteProperty<Derived,Base>::equals(const Property &other) const
{
    return ConcreteProperty<Derived,Base>::operator==(other);
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
void ConcreteProperty<Derived,Base>::save(QDataStream &ds) const
{
    ds << static_cast<const Derived&>(*this);
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
void ConcreteProperty<Derived,Base>::load(QDataStream &ds)
{
    ds >> static_cast<Derived&>(*this);
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
ConcreteProperty<Derived,Base>&
ConcreteProperty<Derived,Base>::operator=(const ConcreteProperty<Derived,Base> &other)
{
    Base::operator=(other);
    return *this;
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
bool ConcreteProperty<Derived,Base>::operator==(
                                    const ConcreteProperty<Derived,Base> &other) const
{
    return Base::operator==(other);
}

template<class Derived, class Base>
SIRE_OUTOFLINE_TEMPLATE
bool ConcreteProperty<Derived,Base>::operator!=(
                                    const ConcreteProperty<Derived,Base> &other) const
{
    return Base::operator!=(other);
}

///////
/////// Implementation of PropPtr<T>
///////

/** Assert that this is sane - this is to make sure that the 
    Property really is derived from 'T' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PropPtr<T>::assertSane() const
{
    const Property &p = PropPtrBase::read();

    if (not p.isA<T>())
    {
        PropPtrBase::throwCastingError( p.what(), T::typeName() );
    }
}

/** Return the null object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T> PropPtr<T>::null()
{
    return PropPtr<T>( T::null() );
}

/** Return whether this is the null object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PropPtr<T>::isNull() const
{
    return PropPtrBase::operator==( T::null() );
}

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::PropPtr() : PropPtrBase( PropPtr<T>::null() )
{}

/** Construct from the object 'obj' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::PropPtr(const T &obj) : PropPtrBase(obj)
{}

/** Construct from a pointer to the passed object 'obj' - this
    takes over ownership of the pointer and can delete it at any time */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::PropPtr(T *obj) : PropPtrBase(obj)
{}

/** Construct from the passed property

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::PropPtr(const Property &property) : PropPtrBase(property)
{
    this->assertSane();
}

/** Construct from a pointer of another type */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::PropPtr(const PropPtrBase &other) : PropPtrBase(other)
{
    this->assertSane();
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::PropPtr(const PropPtr<T> &other) : PropPtrBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::~PropPtr()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>& PropPtr<T>::operator=(const PropPtr<T> &other)
{
    PropPtrBase::operator=(other);
    return *this;
}

/** Create a pointer that points to the object 'obj' (or a copy, if this
    object is on the stack) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>& PropPtr<T>::operator=(const T &obj)
{
    return this->operator=( PropPtr<T>(obj) );
}

/** Create a pointer that points to the object 'obj' - this takes
    over ownership of 'obj' and can delete it at any time */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>& PropPtr<T>::operator=(T *obj)
{
    return this->operator=( PropPtr<T>(obj) );
}

/** Construct from the passed property

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>& PropPtr<T>::operator=(const Property &property)
{
    return this->operator=( PropPtr<T>(property) );
}

/** Construct from the passed pointer to a property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>& PropPtr<T>::operator=(const PropPtrBase &property)
{
    return this->operator=( PropPtr<T>(property) );
}

/** Return a reference to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PropPtr<T>::read() const
{
    //This is only safe as everything in this class ensures that 
    //the base pointer really is a pointer to the derived type 'T' */
    return static_cast<const T&>( PropPtrBase::read() );
}

/** Return an editable reference to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& PropPtr<T>::edit()
{
    //This is only safe as everything in this class ensures that 
    //the base pointer really is a pointer to the derived type 'T' */
    return static_cast<T&>( PropPtrBase::edit() );
}

/** Return the raw pointer to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PropPtr<T>::operator->() const
{
    return &(this->read());
}

/** Return a reference to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PropPtr<T>::operator*() const
{
    return this->read();
}

/** Return the raw pointer to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PropPtr<T>::data() const
{
    return &(this->read());
}

/** Return the raw pointer to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PropPtr<T>::constData() const
{
    return &(this->read());
}

/** Return a modifiable pointer to the data */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PropPtr<T>::data()
{
    return &(this->edit());
}

/** Allow automatic casting to a T() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PropPtr<T>::operator const T&() const
{
    return this->read();
}

///////
/////// Implementation of GlobalPropPtr<T>
///////

/** Assert that this is sane - this is to make sure that the 
    Property really is derived from 'T' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void GlobalPropPtr<T>::assertSane() const
{
    const Property &p = GlobalPropPtrBase::read();

    if (not p.isA<T>())
    {
        GlobalPropPtrBase::throwCastingError( p.what(), T::typeName() );
    }
}

/** Return the null object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T> GlobalPropPtr<T>::null()
{
    return GlobalPropPtr<T>( T::null() );
}

/** Return whether this is the null object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool GlobalPropPtr<T>::isNull() const
{
    return GlobalPropPtrBase::operator==( T::null() );
}

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::GlobalPropPtr() : GlobalPropPtrBase( GlobalPropPtr<T>::null() )
{}

/** Construct from the object 'obj' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::GlobalPropPtr(const T &obj) : GlobalPropPtrBase(obj)
{}

/** Construct from a pointer to the passed object 'obj' - this
    takes over ownership of the pointer and can delete it at any time */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::GlobalPropPtr(T *obj) : GlobalPropPtrBase(obj)
{}

/** Construct from the passed property

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::GlobalPropPtr(const Property &property) : GlobalPropPtrBase(property)
{
    this->assertSane();
}

/** Construct from a pointer of another type */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::GlobalPropPtr(const GlobalPropPtrBase &other) : GlobalPropPtrBase(other)
{
    this->assertSane();
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::GlobalPropPtr(const GlobalPropPtr<T> &other) : GlobalPropPtrBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::~GlobalPropPtr()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>& GlobalPropPtr<T>::operator=(const GlobalPropPtr<T> &other)
{
    GlobalPropPtrBase::operator=(other);
    return *this;
}

/** Create a pointer that points to the object 'obj' (or a copy, if this
    object is on the stack) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>& GlobalPropPtr<T>::operator=(const T &obj)
{
    return this->operator=( GlobalPropPtr<T>(obj) );
}

/** Create a pointer that points to the object 'obj' - this takes
    over ownership of 'obj' and can delete it at any time */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>& GlobalPropPtr<T>::operator=(T *obj)
{
    return this->operator=( GlobalPropPtr<T>(obj) );
}

/** Construct from the passed property

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>& GlobalPropPtr<T>::operator=(const Property &property)
{
    return this->operator=( GlobalPropPtr<T>(property) );
}

/** Construct from the passed pointer to a property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>& GlobalPropPtr<T>::operator=(const GlobalPropPtrBase &property)
{
    return this->operator=( GlobalPropPtr<T>(property) );
}

/** Return a reference to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& GlobalPropPtr<T>::read() const
{
    //This is only safe as everything in this class ensures that 
    //the base pointer really is a pointer to the derived type 'T' */
    return static_cast<const T&>( GlobalPropPtrBase::read() );
}

/** Return the raw pointer to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* GlobalPropPtr<T>::operator->() const
{
    return &(this->read());
}

/** Return a reference to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& GlobalPropPtr<T>::operator*() const
{
    return this->read();
}

/** Return the raw pointer to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* GlobalPropPtr<T>::data() const
{
    return &(this->read());
}

/** Return the raw pointer to the object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* GlobalPropPtr<T>::constData() const
{
    return &(this->read());
}

/** Allow automatic casting to a T() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
GlobalPropPtr<T>::operator const T&() const
{
    return this->read();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

/** This is the specialised pointer that is used to hold a basic
    Property object
    
    @author Christopher Woods
*/
template<>
class SIREBASE_EXPORT PropPtr<Property> : public PropPtrBase
{
public:
    PropPtr();

    PropPtr(const Property &property);
    PropPtr(Property *property);

    PropPtr(const PropPtrBase &other);

    PropPtr(const PropPtr<Property> &other);
    
    ~PropPtr();
    
    PropPtr<Property>& operator=(const Property &property);
    PropPtr<Property>& operator=(Property *property);

    PropPtr<Property>& operator=(const PropPtr<Property> &other);
    PropPtr<Property>& operator=(const PropPtrBase &property);

    const Property* operator->() const;
    const Property& operator*() const;
    
    const Property& read() const;
    Property& edit();

    const Property* data() const;
    const Property* constData() const;

    Property* data();
    
    operator const Property&() const;

    bool isNull() const;

    static PropPtr<Property> null();

    void assertSane() const
    {}
};

/** This is the specialised global pointer that is used to hold a basic
    Property object
    
    @author Christopher Woods
*/
template<>
class SIREBASE_EXPORT GlobalPropPtr<Property> : public GlobalPropPtrBase
{
public:
    GlobalPropPtr();

    GlobalPropPtr(const Property &property);
    GlobalPropPtr(Property *property);

    GlobalPropPtr(const GlobalPropPtrBase &other);

    GlobalPropPtr(const GlobalPropPtr<Property> &other);
    
    ~GlobalPropPtr();
    
    GlobalPropPtr<Property>& operator=(const Property &property);
    GlobalPropPtr<Property>& operator=(Property *property);

    GlobalPropPtr<Property>& operator=(const GlobalPropPtr<Property> &other);
    GlobalPropPtr<Property>& operator=(const GlobalPropPtrBase &property);

    const Property* operator->() const;
    const Property& operator*() const;
    
    const Property& read() const;

    const Property* data() const;
    const Property* constData() const;
    
    operator const Property&() const;

    bool isNull() const;

    static GlobalPropPtr<Property> null();

    void assertSane() const
    {}
};

/** Create the type to hold a pointer to a generic Property */
typedef PropPtr<Property> PropertyPtr;

/** Create the type to hold a global pointer to a generic Property */
typedef GlobalPropPtr<Property> GlobalPropertyPtr;

} // end of namespace SireBase

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Serialise a property pointer to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireBase::PropPtr<T> &prop)
{
    ds << static_cast<const SireBase::PropPtrBase&>(prop);
    return ds;
}

/** Read a property pointer from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireBase::PropPtr<T> &prop)
{
    SireBase::PropPtr<T> ptr;
    
    ds >> static_cast<SireBase::PropPtrBase&>(ptr);
    
    ptr.assertSane();
    
    prop = ptr;
    
    return ds;
}

/** Serialise a property pointer to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireBase::GlobalPropPtr<T> &prop)
{
    ds << static_cast<const SireBase::GlobalPropPtrBase&>(prop);
    return ds;
}

/** Read a property pointer from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireBase::GlobalPropPtr<T> &prop)
{
    SireBase::GlobalPropPtr<T> ptr;
    
    ds >> static_cast<SireBase::GlobalPropPtrBase&>(ptr);
    
    ptr.assertSane();
    
    prop = ptr;
    
    return ds;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireBase::NullProperty);
Q_DECLARE_METATYPE(SireBase::VariantProperty);

SIRE_EXPOSE_CLASS( SireBase::Property )
SIRE_EXPOSE_CLASS( SireBase::NullProperty )
SIRE_EXPOSE_CLASS( SireBase::VariantProperty )

SIRE_EXPOSE_PROPERTY( SireBase::PropertyPtr, SireBase::Property )

SIRE_END_HEADER

#endif
