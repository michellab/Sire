/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREBASE_PACKEDARRAY2D_HPP
#define SIREBASE_PACKEDARRAY2D_HPP

#include "qvariant_metatype.h"

#include <QVarLengthArray>

#include "packedarray2d.h"
#include "quickcopy.hpp"

#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <boost/assert.hpp>

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class PackedArray2D;

namespace detail
{
template<class T>
class PackedArray2D_Array;
}

}

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::PackedArray2D<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::PackedArray2D<T>&);

template<class T>
QDataStream& operator<<(QDataStream&, 
                        const SireBase::detail::PackedArray2D_Array<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, 
                        SireBase::detail::PackedArray2D_Array<T>&);

namespace SireBase
{

template<class T>
class PackedArray2D;

namespace detail
{

template<class T>
class SharedArray2DPtr;

template<class T>
class PackedArray2DData;

template<class T>
class PackedArray2D_ArrayData;

template<class T>
class PackedArray2D_Array;

/** This is a simple class that is used to manage the memory and 
    reference counting for the PackedArray2D */
template<class T>
class SIREBASE_EXPORT PackedArray2DMemory : public PackedArray2DMemoryBase
{
public:
    static char* create(quint32 narrays, quint32 nvalues);
    static char* detach(char *this_ptr, quint32 this_idx);
    static void destroy(PackedArray2DData<T> *arrays);
    
    static void incref(char *this_ptr, quint32 this_idx);
    static void decref(char *this_ptr, quint32 this_idx);

    static SharedArray2DPtr< PackedArray2DData<T> > shared_null;
};

/** This is the implicitly shared pointer class that 
    is used to hold any of the PackedArray2DMemory allocated objects
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT SharedArray2DPtr
{

public:
    SharedArray2DPtr() : ptr(0)
    {}

    SharedArray2DPtr(const T *p)
    {
        ptr = const_cast<T*>(p);
    
        if (ptr)
            ptr->incref();
    }
    
    SharedArray2DPtr(const SharedArray2DPtr &other) : ptr(other.ptr)
    {
        if (ptr)
            ptr->incref();
    }
    
    SharedArray2DPtr(SharedArray2DPtr &&other) : ptr(other.ptr)
    {
        other.ptr = 0;
    }
    
    ~SharedArray2DPtr()
    {
        if (ptr)
            ptr->decref();
    }
    
    SharedArray2DPtr<T>& operator=(const SharedArray2DPtr &other)
    {
        if (ptr != other.ptr)
        {
            T *new_ptr = other.ptr;
            
            //increment the other reference count
            if (new_ptr)
                new_ptr->incref();
                
            //decrement our reference count
            if (ptr)
                ptr->decref();
                
            //set the new pointer
            ptr = new_ptr;
        }
        
        return *this;
    }
    
    SharedArray2DPtr<T>& operator=(SharedArray2DPtr<T> &&other)
    {
        T *new_ptr = other.ptr;
        other.ptr = 0;
        
        if (ptr)
            ptr->decref();
        
        ptr = new_ptr;
        
        return *this;
    }
    
    const T& operator*() const
    {
        return *ptr;
    }
    
    const T* operator->() const
    {
        return ptr;
    }
    
    T& operator*()
    {
        if (ptr)
            ptr = ptr->detach();
            
        return *ptr;
    }
    
    T* operator->()
    {
        if (ptr)
            ptr = ptr->detach();
            
        return ptr;
    }
    
    const T* data() const
    {
        return ptr;
    }
    
    const T* constData() const
    {
        return ptr;
    }
    
    T* data()
    {
        if (ptr)
            ptr = ptr->detach();
        
        return ptr;
    }
    
    /** Detach this pointer */
    void detach()
    {
        if (ptr)
            ptr = ptr->detach();
    }
    
    /** Assign this pointer to point at 'weakptr' 
        but *without* changing the reference count.
        You ABSOLUTELY MUST ensure that you call 
        SharedArray2DPtr::weakRelease() before this 
        pointer is deleted or reassigned, so 
        as to not decrement the reference count incorrectly! */
    void weakAssign(T *weakptr)
    {
        if (ptr)
            ptr->decref();
            
        ptr = weakptr;
    }
    
    /** Release the pointer *without* decrementing the
        reference count. You should only call this
        function if the pointer was assigned using 
        the 'weakAssign()' function */
    void weakRelease()
    {
        ptr = 0;
    }
    
private:
    /** Actual pointer */
    T *ptr;
};

/** This class is used to hold all of the metadata about the 
    packed array of arrays
*/
template<class T>
class SIREBASE_EXPORT PackedArray2DData : public PackedArray2DDataBase
{

friend class PackedArray2DMemory<T>;

public:
    PackedArray2DData();
    PackedArray2DData(quint32 narrays, quint32 nvalues);
    
    PackedArray2DData(const PackedArray2DData &other);
    
    ~PackedArray2DData();

    void decref();
    
    PackedArray2DData* detach();
    
    const PackedArray2D_ArrayData<T>* nullArray() const;
    
    const PackedArray2D_ArrayData<T>* arrayDataData() const;
    const PackedArray2D_Array<T>* arrayData() const;

    const T* valueData() const;

    PackedArray2D_ArrayData<T>* arrayDataData();
    PackedArray2D_Array<T>* arrayData();

    T* valueData();
    
    void setNValuesInArray(quint32 i, quint32 nvalues);
    
    void close();
};

/** This class holds the metadata about an individual array
    in the packed collection of arrays */
template<class T>
class SIREBASE_EXPORT PackedArray2D_ArrayData : public PackedArray2D_ArrayDataBase
{

friend class PackedArray2D<T>;
friend class PackedArray2DData<T>;
friend class PackedArray2DMemory<T>;

public:
    PackedArray2D_ArrayData();
    PackedArray2D_ArrayData(quint32 this_idx);
    
    PackedArray2D_ArrayData(const PackedArray2D_ArrayData &other);
    
    ~PackedArray2D_ArrayData();
    
    void incref();
    void decref();
    
    PackedArray2D_ArrayData<T>* detach();
    
    PackedArray2DData<T>* extract() const;
    
    const T* valueData() const;
    T* valueData();
};

/** This is the implicitly shared array class that is packed
    together to make a PackedArray2D<T>. This class is mainly used
    via its alias - PackedArray2D<T>::Array. Note that the size
    of a PackedArray2D<T>::Array is fixed at construction. You
    cannot change the size as this array is packed contiguously
    in memory with other arrays.
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT PackedArray2D_Array
{

friend class PackedArray2DMemory<T>;
friend class PackedArray2D<T>;

friend QDataStream& ::operator<<<>(QDataStream&, const PackedArray2D_Array<T>&);
friend QDataStream& ::operator>><>(QDataStream&, PackedArray2D_Array<T>&);

public:
    PackedArray2D_Array();
    
    PackedArray2D_Array(quint32 sz);
    PackedArray2D_Array(quint32 sz, const T &value);
    
    PackedArray2D_Array(const QVector<T> &values);
    
    PackedArray2D_Array(const PackedArray2D_Array<T> &other);
    
    ~PackedArray2D_Array();
    
    PackedArray2D_Array<T>& operator=(const PackedArray2D_Array<T> &other);
    
    bool operator==(const PackedArray2D_Array<T> &other) const;
    bool operator!=(const PackedArray2D_Array<T> &other) const;

    const T& operator[](quint32 i) const;
    T& operator[](quint32 i);
    
    const T& at(quint32 i) const;

    QString toString() const;

    int count() const;
    int size() const;
    
    int nValues() const;
    
    bool isEmpty() const;
    
    const T* data() const;
    T* data();
    
    const T* constData() const;
    
    void update(const PackedArray2D_Array<T> &other);
    
    QVector<T> toQVector() const;
    
    void assertValidIndex(quint32 i) const;

protected:
    PackedArray2D_Array(detail::PackedArray2D_ArrayData<T> *data);
    
private:
    /** Implicitly shared pointer to the array */
    detail::SharedArray2DPtr< detail::PackedArray2D_ArrayData<T> > d;
};

}

/** This class provides an array of arrays of type T, where
    each array is packed in memory so that each object of type
    T is contiguous. This allows rapid indexing over each
    object in each array, or over all objects in all arrays.
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT PackedArray2D
{

friend class detail::PackedArray2DMemory<T>;

friend QDataStream& ::operator<<<>(QDataStream&, const PackedArray2D<T>&);
friend QDataStream& ::operator>><>(QDataStream&, PackedArray2D<T>&);

public:
    typedef typename detail::PackedArray2D_Array<T> Array;
    typedef T value_type;

    PackedArray2D();
    
    PackedArray2D(const Array &array);
    PackedArray2D(const QVector<Array> &arrays);
    
    PackedArray2D(const QVector<T> &values);
    PackedArray2D(const QVector< QVector<T> > &values);
    
    PackedArray2D(const PackedArray2D<T> &array0, const PackedArray2D<T> &array1);
    
    PackedArray2D(const PackedArray2D<T> &other);
    
    ~PackedArray2D();
    
    PackedArray2D<T>& operator=(const PackedArray2D<T> &other);
    
    bool operator==(const PackedArray2D<T> &other) const;
    bool operator!=(const PackedArray2D<T> &other) const;
    
    const Array& operator[](quint32 i) const;
    
    const T& operator()(quint32 i, quint32 j) const;
    T& operator()(quint32 i, quint32 j);
    
    const Array& at(quint32 i) const;
    const T& at(quint32 i, quint32 j) const;
    
    int count() const;
    int size() const;
    
    int nArrays() const;
    int nValues() const;
    int nValues(quint32 i) const;

    bool isEmpty() const;
    
    QString toString() const;

    void detach();

    const Array* data() const;
    const Array* constData() const;
    
    const T* data(quint32 i) const;
    T* data(quint32 i);
    
    const T* constData(quint32 i) const;
    
    const T* valueData() const;
    T* valueData();
    
    const T* constValueData() const;

    QVector<T> toQVector() const;
    QVector< QVector<T> > toQVectorVector() const;

    void update(quint32 i, const Array &array);
    void update(quint32 i, const QVector<T> &array);

    template<class C>
    void updateAll(const C &idxs, const PackedArray2D<T> &arrays);
    
    void updateAll(const QVarLengthArray<int> &idxs,    
                   const PackedArray2D<T> &arrays);
    
    template<class C>
    void updateAll(const C &idxs, const QVector< QVector<T> > &arrays);
    
    void append(const Array &array);
    void append(const PackedArray2D<T> &arrays);
    
    void append(const QVector<T> &array);
    void append(const QVector< QVector<T> > &arrays);
    
    void remove(quint32 i);
    
    template<class C>
    void removeAll(const C &idxs);

    void removeAll(const QVarLengthArray<int> &idxs);

    void assertValidIndex(quint32 i) const;

    PackedArray2D<QVariant> toVariant() const;
    
    static PackedArray2D<T> fromVariant(const PackedArray2D<QVariant> &variant);

private:
    /** Implicitly shared pointer to the array data */
    detail::SharedArray2DPtr< detail::PackedArray2DData<T> > d;
};

}

namespace SireBase
{

namespace detail
{

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

////////
//////// Implementation of PackedArray2DMemory
////////

/** Pointer to the global null PackedArray */
template<class T>
SharedArray2DPtr< PackedArray2DData<T> > PackedArray2DMemory<T>::shared_null;

/** Create space for the arrays. The layout in memory is;

    --------------------------------------------------
    | A |   B   |    C    |            D             |
    --------------------------------------------------

      A = PackedArray2DData<T> object          1 * sizeof(PackedArray2DData<T>)
      B = PackedArray2D_Array<T> objects      narrays * sizeof(PackedArray2D_Array<T>)
      C = PackedArray2D_ArrayData<T> objects   narrays * sizeof( " _ArrayData<T>)
      D = Array of T objects                 nvalues * sizeof(T)
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
char* PackedArray2DMemory<T>::create(quint32 narrays, quint32 nvalues)
{
    //create space for the data
    char *storage = PackedArray2DMemoryBase::create(narrays, nvalues, 
                                            sizeof(PackedArray2DData<T>),
                                            sizeof(PackedArray2D_Array<T>),
                                            sizeof(PackedArray2D_ArrayData<T>),
                                            sizeof(T));
    
    try
    {
        quint32 sz = PackedArray2DMemoryBase::getSize(narrays, nvalues,   
                                            sizeof(PackedArray2DData<T>),
                                            sizeof(PackedArray2D_Array<T>),
                                            sizeof(PackedArray2D_ArrayData<T>),
                                            sizeof(T));
    
        //now we need to construct each object in turn at its correct
        //location in this space
        PackedArray2DData<T> *arraydata 
                    = new (storage) PackedArray2DData<T>(narrays, nvalues);
                    
        //advance the index into the storage array so that we are now
        //pointing just after the PackedArray2DData
        quint32 idx = sizeof(PackedArray2DData<T>);
                    
        //this is the location of the first PackedArray2D_Array<T> :-)
        PackedArray2DMemoryBase::setArray0(arraydata, idx);
        quint32 array0_idx = idx;
        
        //the first PackedArray2D_ArrayData lies at this 
        //index + narrays*sizeof(PackedArray2D_Array<T>)
        quint32 dataidx = idx + narrays * sizeof(PackedArray2D_Array<T>);
        
        PackedArray2DMemoryBase::setArrayData0(arraydata, dataidx);

        if (narrays > 0)
        {
            //loop over each array and create it in its place
            for (quint32 i=0; i<narrays; ++i)
            {   
                //assert that there is sufficient space in the array
                BOOST_ASSERT(idx + sizeof(PackedArray2D_Array<T>) < sz);
                BOOST_ASSERT(dataidx + sizeof(PackedArray2D_ArrayData<T>) < sz);
            
                //create the ArrayData, letting it know where it
                //is relative to the beginning of the storage array
                PackedArray2D_ArrayData<T> *data 
                        = new (storage + dataidx) PackedArray2D_ArrayData<T>(dataidx);
                
                //now create the PackedArray2D_Array<T> that uses this data
                new (storage + idx) PackedArray2D_Array<T>(data);
        
                //advance the index into the storage array to point
                //just after the just-created PackedArray2D_ArrayData
                idx += sizeof(PackedArray2D_Array<T>);
                dataidx += sizeof(PackedArray2D_ArrayData<T>);
            }
            
            //make idx point just after all of the PackedArray2D_ArrayDatas...
            BOOST_ASSERT( idx == array0_idx 
                                + narrays*sizeof(PackedArray2D_Array<T>) );
                                
            idx = dataidx;
        }
        else
        {
            //we need to create space for the null PackedArray2D_Array<T>
            quint32 data_idx = idx + sizeof(PackedArray2D_Array<T>);
            PackedArray2D_ArrayData<T> *array_arraydata 
                        = new (storage + data_idx) PackedArray2D_ArrayData<T>(data_idx);
            
            new (storage + idx) PackedArray2D_Array<T>(array_arraydata);
            
            PackedArray2DMemoryBase::setNValues(array_arraydata, 0);
            PackedArray2DMemoryBase::setValue0(array_arraydata, 0);
            
            idx += sizeof(PackedArray2D_Array<T>) 
                    + sizeof(PackedArray2D_ArrayData<T>);
        }
        
        //we are now at the location of the first item in the array
        PackedArray2DMemoryBase::setValue0(arraydata, idx);
        
        //only call constructors if this is a complex type
        if (QTypeInfo<T>::isComplex)
        {
            //loop over each object and create it in place
            for (quint32 i=0; i<nvalues; ++i)
            {
                BOOST_ASSERT(idx + sizeof(T) <= sz);

                new (storage + idx) T();
                    
                idx += sizeof(T);
            }
        }
        else
            idx += nvalues * sizeof(T);
        
        //we should now be at the end of the storage
        BOOST_ASSERT( idx == sz );
                
        return storage;
    }
    catch(...)
    {
        delete[] storage;
        throw;
    }
}

/** Increase the reference count */
template<class T>
SIRE_INLINE_TEMPLATE
void PackedArray2DMemory<T>::incref(char *this_ptr, quint32 this_idx)
{
    ( (PackedArray2DData<T>*)(PackedArray2DMemoryBase::getRoot(this_ptr,this_idx)) )
                    ->incref();
}

/** Decrease the reference count */
template<class T>
SIRE_INLINE_TEMPLATE
void PackedArray2DMemory<T>::decref(char *this_ptr, quint32 this_idx)
{
    ( (PackedArray2DData<T>*)(PackedArray2DMemoryBase::getRoot(this_ptr,this_idx)) )
                    ->decref();
}

/** Destroy the object 'array' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DMemory<T>::destroy(PackedArray2DData<T> *array)
{
    //we need to delete it in the opposite order to creation - first
    //lets delete all of the objects
    quint32 nvalues = array->nValues();

    char *storage = (char*)array;
        
    //only call destructors if this is a complex type
    if (nvalues > 0 and QTypeInfo<T>::isComplex)
    {
        T *values = array->valueData();
        
        for (qint32 i=nvalues-1; i>=0; --i)
        {
            values[i].~T();
        }
    }
        
    //now delete all of the arrays
    quint32 narrays = array->nArrays();
        
    if (narrays > 0)
    {
        PackedArray2D_Array<T> *arrays = array->arrayData();
            
        for (qint32 i=narrays-1; i>=0; --i)
        {
            PackedArray2D_ArrayData<T> *arraydata 
                 = const_cast<PackedArray2D_ArrayData<T>*>( arrays[i].d.constData() );
            
            //delete the PackedArray2D_ArrayData
            arraydata->~PackedArray2D_ArrayData<T>();
            
            //remove the shared pointer
            arrays[i].d.weakRelease();
            
            //delete the PackedArray2D_Array<T>
            arrays[i].~PackedArray2D_Array<T>();
        }
    }
        
    //delete the PackedArray2DData object itself
    array->~PackedArray2DData<T>();
        
    //finally release the memory
    delete[] storage;
}

/** Detach the data pointed to by 'this_ptr' from shared storage */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
char* PackedArray2DMemory<T>::detach(char *this_ptr, quint32 this_idx)
{
    //get a pointer to the start of the storage for this container
    char *storage = getRoot(this_ptr, this_idx);
    
    //The PackedArray2DData object is at the beginning of this storage array
    PackedArray2DData<T> *arraydata = (PackedArray2DData<T>*) storage;
    
    if (not arraydata->ref.testAndSetRelaxed(1,1))
    {
        //there is more than one reference to this data - it will have to 
        //be cloned - get the size of memory to be cloned
        int sz = getSize(arraydata->nArrays(), arraydata->nValues(),
                         sizeof(PackedArray2DData<T>), 
                         sizeof(PackedArray2D_Array<T>),
                         sizeof(PackedArray2D_ArrayData<T>),
                         sizeof(T));
                         
        //create space for the copy
        char *new_storage = new char[sz];
        
        //copy the data
        quickCopy<char>( new_storage, storage, sz );

        //the first part of the data is the PackedArray2DData object
        PackedArray2DData<T> *new_arraydata = (PackedArray2DData<T>*) new_storage;

        //call the copy constructor if this is a complex type
        if (QTypeInfo<T>::isComplex)
        {
            T *new_array = new_arraydata->valueData();
            const T *old_array = arraydata->valueData();
            
            int nvals = arraydata->nValues();
            
            for (int i=0; i<nvals; ++i)
            {
                new (new_array) T(*old_array);
                
                ++new_array;
                ++old_array;
            }
        }
        
        //set the reference count of this copy to 1
        new_arraydata->ref = QAtomicInt(1);

        //now loose a reference to the original
        PackedArray2DMemory<T>::decref(this_ptr, this_idx);
        
        //the final step is to update all of the PackedArray2D_Array<T>
        //pointers that exist in this array (otherwise they will all
        //point to the original!)
        
        if (new_arraydata->nArrays() > 0)
        {
            PackedArray2D_Array<T> *arrays = new_arraydata->arrayData();
            PackedArray2D_ArrayData<T> *arrays_data = new_arraydata->arrayDataData();
        
            for (quint32 i=0; i < new_arraydata->nArrays(); ++i)
            {
                arrays[i].d.weakRelease();
                arrays[i].d.weakAssign( &(arrays_data[i]) );
            }
        }
        else
        {
            //update the null PackedArray2D_Array<T>
            PackedArray2D_Array<T> *arrays = new_arraydata->arrayData();
            PackedArray2D_ArrayData<T> *arrays_data = new_arraydata->arrayDataData();

            arrays[0].d.weakRelease();
            arrays[0].d.weakAssign( &(arrays_data[0]) );
        }

        //return a pointer to the clone
        return new_storage + this_idx;
    }
    else
    {
        //only one reference, so no need to clone
        return this_ptr;
    }
}

///////
/////// Implementation of PackedArray2DData<T>
///////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::PackedArray2DData() : PackedArray2DDataBase()
{}

/** Construct to hold the specified number of arrays and values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::PackedArray2DData(quint32 narrays, quint32 nvalues)
                   : PackedArray2DDataBase(narrays, nvalues)
{}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::PackedArray2DData(const PackedArray2DData &other)
                   : PackedArray2DDataBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::~PackedArray2DData()
{}

/** Decrease the reference count of this object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DData<T>::decref()
{
    if (not this->ref.deref())
        PackedArray2DMemory<T>::destroy(this);
}

/** Return a pointer to a copy of this array that has been 
    detached from shared storage */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>* PackedArray2DData<T>::detach()
{
    return (PackedArray2DData<T>*)( PackedArray2DMemory<T>::detach( (char*)(this), 0 ) );
}

/** Return a pointer to the null array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D_ArrayData<T>* PackedArray2DData<T>::nullArray() const
{
    BOOST_ASSERT( this->nArrays() == 0 );
    BOOST_ASSERT( this->getArrayData0() != 0 );
    
    return (const PackedArray2D_ArrayData<T>*)
                    ( this->memory() + this->getArrayData0() );
}

/** Return a raw pointer to the array of PackedArray2D_ArrayData objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D_ArrayData<T>* PackedArray2DData<T>::arrayDataData() const
{
    if (this->getArrayData0() == 0)
        return 0;
    else
        return (const PackedArray2D_ArrayData<T>*)
                    ( this->memory() + this->getArrayData0() );
}

/** Return a raw pointer to the array of PackedArray2D_Array<T> objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D_Array<T>* PackedArray2DData<T>::arrayData() const
{
    if (this->getArray0() == 0)
        return 0;
    else
        return (const PackedArray2D_Array<T>*)( this->memory() + this->getArray0() );
}

/** Return a raw pointer to the array of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2DData<T>::valueData() const
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (const T*)( this->memory() + this->getValue0() );
}

/** Return a raw pointer to the array of PackedArray2D_ArrayData objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>* PackedArray2DData<T>::arrayDataData()
{
    if (this->getArrayData0() == 0)
        return 0;
    else
        return (PackedArray2D_ArrayData<T>*)( this->memory() + this->getArrayData0() );
}

/** Return a raw pointer to the array of PackedArray2D_Array<T> objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>* PackedArray2DData<T>::arrayData()
{
    if (this->getArray0() == 0)
        return 0;
    else
        return (PackedArray2D_Array<T>*)( this->memory() + this->getArray0() );
}

/** Return a raw pointer to the array of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2DData<T>::valueData()
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (T*)( this->memory() + this->getValue0() );
}

/** Set the number of values in the ith array - can only be called
    before the close() (so don't call it after closing the arrays) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DData<T>::setNValuesInArray(quint32 i, quint32 nvalues)
{
    BOOST_ASSERT( i < this->nArrays() );
    PackedArray2DMemoryBase::setNValues(&(this->arrayDataData()[i]), nvalues);
}

/** Close out this array - this fixes the number of values
    in each array in this array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DData<T>::close()
{
    PackedArray2D_ArrayData<T> *arraydata = this->arrayDataData();
    
    //ensure that everything adds up
    quint32 n_assigned_vals = 0;
    
    quint32 value0 = this->getValue0();
    
    //close out each array
    for (quint32 i=0; i<this->nArrays(); ++i)
    {
        PackedArray2D_ArrayData<T> &array = arraydata[i];
        
        //tell the array where its first value is located
        if (array.nValues() > 0)
            PackedArray2DMemoryBase::setValue0(&array, 
                                        value0 + n_assigned_vals * sizeof(T));
        else
            PackedArray2DMemoryBase::setValue0(&array, 0);

        n_assigned_vals += array.nValues();
    }

    if (n_assigned_vals != this->nValues())
        throw SireError::program_bug( QObject::tr(
                "Serious program bug - n_assigned_vals (%1) does not equal "
                "this->nValues() (%2).")
                    .arg(n_assigned_vals).arg(this->nValues()), CODELOC );
}

////////
//////// Implementation of PackedArray2D_ArrayData
////////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::PackedArray2D_ArrayData()
                         : PackedArray2D_ArrayDataBase()
{}

/** Construct the data, telling it that it is at 
    index 'this_idx' in the memory storage array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::PackedArray2D_ArrayData(quint32 this_idx)
                         : PackedArray2D_ArrayDataBase(this_idx)
{}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::PackedArray2D_ArrayData(const PackedArray2D_ArrayData<T> &other)
                         : PackedArray2D_ArrayDataBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::~PackedArray2D_ArrayData()
{}

/** Return a pointer to a PackedArray2D_ArrayData that has been detached from 
    shared storage (could be a pointer to this PackedArray2D_ArrayData) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>* PackedArray2D_ArrayData<T>::detach()
{
    return (PackedArray2D_ArrayData<T>*)
           ( PackedArray2DMemory<T>::detach( (char*)this, this->getThisArray() ) );
}

/** Return a pointer to a PackedArray2DData object that contains 
    just the data of this PackedArray2D_Array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>* PackedArray2D_ArrayData<T>::extract() const
{
    const PackedArray2DData<T> *array = (const PackedArray2DData<T>*)(this->memory());
    
    if (array->nArrays() <= 1)
        //this is already extracted!
        return const_cast<PackedArray2DData<T>*>(array);
    
    char *new_storage = PackedArray2DMemory<T>::create(1, this->nValues());
    
    try
    {
                
    //ok, we need to extract!
    PackedArray2DData<T> *new_array = (PackedArray2DData<T>*)(new_storage);

    new_array->setNValuesInArray(0, this->nValues());
    new_array->close();
    
    //copy the objects
    quickCopy<T>(new_array->valueData(),
                 this->valueData(), this->nValues());

    return new_array;

    }
    catch(...)
    {
        delete[] new_storage;
    }
    
    return 0;
}

/** Increment the reference count for this object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_ArrayData<T>::incref()
{
    PackedArray2DMemory<T>::incref( (char*)this, this->getThisArray() );
}

/** Decrement the reference count for this object - this may delete it! */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_ArrayData<T>::decref()
{
    PackedArray2DMemory<T>::decref( (char*)this, this->getThisArray() );
}

/** Return a pointer to the first object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D_ArrayData<T>::valueData() const
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (const T*)( this->memory() + this->getValue0() );
}

/** Return a pointer to the first object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D_ArrayData<T>::valueData()
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (T*)( this->memory() + this->getValue0() );
}

///////////
/////////// Implementation of PackedArray2D_Array
///////////

template<class T>
static const SharedArray2DPtr< PackedArray2DData<T> >& getSharedNull()
{
    if (detail::PackedArray2DMemory<T>::shared_null.constData() == 0)
    {
        detail::PackedArray2DMemory<T>::shared_null 
                = (PackedArray2DData<T>*)( PackedArray2DMemory<T>::create(0,0) );

        detail::PackedArray2DMemory<T>::shared_null->close();
    }

    return detail::PackedArray2DMemory<T>::shared_null;
}

template<class T>
static SharedArray2DPtr< PackedArray2D_ArrayData<T> > createNullArray()
{
    const SharedArray2DPtr< PackedArray2DData<T> > &array 
                            = SireBase::detail::getSharedNull<T>();
    
    return SharedArray2DPtr< PackedArray2D_ArrayData<T> >( array->nullArray() );
}

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array()
                       : d( SireBase::detail::createNullArray<T>() )
{}

template<class T>
static SharedArray2DPtr< PackedArray2D_ArrayData<T> > createArray(quint32 sz)
{
    if (sz == 0)
        return SireBase::detail::createNullArray<T>();
        
    //construct space for 1 array of sz objects
    char *storage = PackedArray2DMemory<T>::create(1, sz);

    PackedArray2DData<T> *array = (PackedArray2DData<T>*)storage;

    array->setNValuesInArray(0, sz);

    array->close();

    return SharedArray2DPtr< PackedArray2D_ArrayData<T> >( &(array->arrayDataData()[0]) );
}

/** Construct an array of size 'sz'. All of the objects will be initialised  
    with their default values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(quint32 sz)
                     : d( SireBase::detail::createArray<T>(sz) )
{}

/** Construct an array of size 'sz' where all of the objects have value 'value'. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(quint32 sz, const T &value)
                     : d( SireBase::detail::createArray<T>(sz) )
{
    if (sz == 0)
        return;

    T *data = d->valueData();

    //initialise all of the values
    for (quint32 i=0; i<sz; ++i)
    {
        data[i] = value;
    }
}

/** Construct from the passed array of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(const QVector<T> &values)
{
    PackedArray2D<T> array(values);
    this->operator=( array.constData()[0] );
}

/** Private constructor used internally */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(detail::PackedArray2D_ArrayData<T> *data)
{
    d.weakAssign(data);
}

/** Copy constructor - fast as this class is implicitly shared */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(const PackedArray2D_Array<T> &other)
                     : d(other.d)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::~PackedArray2D_Array()
{}

/** Copy assignment operator - fast as this class is implicitly shared.
    However, as this copies a view of the whole packed array, you cannot
    use this to assign the contents of this array within its current 
    PackedArray2D - to do this, you must use the 'update()' function. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>& PackedArray2D_Array<T>::operator=(
                                        const PackedArray2D_Array<T> &other)
{
    this->d = other.d;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D_Array<T>::operator==(const PackedArray2D_Array<T> &other) const
{
    if (this->d.constData() == other.d.constData())
        return true;

    if (this->nValues() != other.nValues())
        return false;
        
    const T *this_data = this->constData();
    const T *other_data = other.constData();
    
    for (int i=0; i<this->nValues(); ++i)
    {
        if (this_data[i] != other_data[i])
            return false;
    }
    
    return true;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D_Array<T>::operator!=(const PackedArray2D_Array<T> &other) const
{
    return not this->operator==(other);
}

/** Assert that the index 'i' is valid for this array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_Array<T>::assertValidIndex(quint32 i) const
{
    if (i >= quint32(this->nValues()))
        detail::throwPackedArray2D_Array_invalidIndex(i, this->nValues());
}

/** Return a raw pointer to the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D_Array<T>::data() const
{
    return this->d->valueData();
}

/** Return a modifiable raw pointer to the array
     *DO NOT* delete this pointer! */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D_Array<T>::data()
{
    return this->d->valueData();
}

/** Return a raw pointer to the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D_Array<T>::constData() const
{
    return this->d->valueData();
}

/** Return a reference to the ith object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D_Array<T>::operator[](quint32 i) const
{
    this->assertValidIndex(i);
    return this->data()[i];
}

/** Return a reference to the ith object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& PackedArray2D_Array<T>::operator[](quint32 i)
{
    this->assertValidIndex(i);
    return this->data()[i];
}

/** Return a reference to the ith object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D_Array<T>::at(quint32 i) const
{
    this->assertValidIndex(i);
    return this->data()[i];
}

/** Return the number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D_Array<T>::count() const
{
    return d->nValues();
}

/** Return the number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D_Array<T>::size() const
{
    return this->d->nValues();
}

/** Return the number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D_Array<T>::nValues() const
{
    return this->d->nValues();
}

/** Return whether or not this array is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D_Array<T>::isEmpty() const
{
    return this->count() == 0;
}

/** Return the QVector that has the same contents as this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVector<T> PackedArray2D_Array<T>::toQVector() const
{
    if (this->isEmpty())
        return QVector<T>();
        
    QVector<T> ret( this->count() );
    
    quickCopy<T>(ret.data(), this->data(), this->count());
    
    return ret;
}

/** Return a string representation of this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString PackedArray2D_Array<T>::toString() const
{
    return Sire::toString( PackedArray2D_Array<T>::toQVector() );
}

/** Update this array so that it has the same contents as 'other'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_Array<T>::update(const PackedArray2D_Array<T> &other)
{
    if (this->size() != other.size())
        detail::throwPackedArray2D_Array_incompatibleError(this->size(),other.size());
        
    if (this->size() == 0)
        return;

    quickCopy<T>(this->data(), other.data(), this->size());
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace detail

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/////////
///////// Implementation of PackedArray2D
/////////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D() : d( SireBase::detail::getSharedNull<T>() )
{}

/** Construct from a passed Array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const typename PackedArray2D<T>::Array &array)
               : d( array.d->extract() )
{}

/** Copy assignment operator - fast as this class is implicitly shared */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>& PackedArray2D<T>::operator=(const PackedArray2D<T> &other)
{
    d = other.d;
    return *this;
}

namespace detail
{

template<class T>
static SharedArray2DPtr< PackedArray2DData<T> > 
createArray(quint32 narrays, quint32 nvals)
{
    if (narrays == 0)
        return SireBase::detail::getSharedNull<T>();
        
    //construct space for narrays arrays of nvals objects
    char *storage = PackedArray2DMemory<T>::create(narrays, nvals);
        
    PackedArray2DData<T> *array = (PackedArray2DData<T>*)storage;
    
    return SharedArray2DPtr< PackedArray2DData<T> >(array);
}

}

/** Construct from an array of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const QVector<T> &values)
                 : d( SireBase::detail::getSharedNull<T>() )
{
    if (values.isEmpty())
        return;

    //count the number of values
    quint32 nvals = values.count();
    
    d = SireBase::detail::createArray<T>(1, nvals);
    
    detail::PackedArray2DData<T> *dptr = d.data();
    T *new_values = dptr->valueData();
    
    //dimension the packed array
    dptr->setNValuesInArray(0, nvals);
    dptr->close();
    
    //now copy all of the data
    quickCopy(new_values, values.constData(), nvals);
}

/** Construct from an array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const QVector<typename PackedArray2D<T>::Array> &arrays)
                 : d( SireBase::detail::getSharedNull<T>() )
{
    if (arrays.isEmpty())
        return;

    else if (arrays.count() == 1)
    {
        this->operator=( PackedArray2D<T>(arrays.at(0)) );
        return;
    }

    //count the number of values and groups...
    quint32 nvals = 0;
    
    const Array *arrays_data = arrays.constData();
    quint32 narrays = arrays.count();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        nvals += arrays_data[i].count();
    }
    
    d = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = d.data();
    T *values = dptr->valueData();
    
    //dimension each packed array
    for (quint32 i=0; i<narrays; ++i)
    {
        //dimension the array
        dptr->setNValuesInArray(i, arrays_data[i].count());
    }
    
    dptr->close();
    
    //now copy all of the data
    for (quint32 i=0; i<narrays; ++i)
    {
        const Array &array = arrays_data[i];
        const T *array_values = array.constData();
        
        quickCopy<T>(values, array_values, array.count());
        
        values += array.count();
    }
}

/** Construct from an array of array of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const QVector< QVector<T> > &values)
                 : d( SireBase::detail::getSharedNull<T>() )
{
    if (values.isEmpty())
        return;
        
    else if (values.count() == 1)
    {
        this->operator=( PackedArray2D<T>(values.at(0)) );
        return;
    }

    //count the number of values and groups...
    quint32 nvals = 0;
    
    const QVector<T> *values_data = values.constData();
    quint32 narrays = values.count();

    for (quint32 i=0; i<narrays; ++i)
    {
        nvals += values_data[i].count();
    }
    
    d = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = d.data();
    
    //dimension each packed array
    for (quint32 i=0; i<narrays; ++i)
    {
        //dimension the array
        dptr->setNValuesInArray(i, values_data[i].count());
    }
    
    dptr->close();
    
    //now copy all of the data
    T *values_array = dptr->valueData();

    for (quint32 i=0; i<narrays; ++i)
    {
        const QVector<T> &array = values_data[i];

        const T *array_values = array.constData();
        
        quickCopy(values_array, array_values, array.count());
        
        values_array += array.count();
    }
}

/** Construct by merging together the passed two arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const PackedArray2D<T> &array0,
                                const PackedArray2D<T> &array1)
{
    if (array0.isEmpty())
    {
        d = array1.d;
        return;
    }
    else if (array1.isEmpty())
    {
        d = array0.d;
        return;
    }

    //count the number of values and groups...
    quint32 nvals = array0.nValues() + array1.nValues();
    quint32 narrays = array0.nArrays() + array1.nArrays();
    
    d = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = d.data();

    const typename PackedArray2D<T>::Array *array0_data = array0.constData();
    const typename PackedArray2D<T>::Array *array1_data = array1.constData();
    
    //dimension each packed array
    for (int i=0; i<array0.nArrays(); ++i)
    {
        //dimension the array
        dptr->setNValuesInArray(i, array0_data[i].count());
    }
    
    for (int i=0; i<array1.nArrays(); ++i)
    {
        dptr->setNValuesInArray(array0.nArrays() + i, array1_data[i].count());
    }
    
    dptr->close();
    
    //now copy all of the data
    T *values_array = dptr->valueData();

    quickCopy(values_array, array0.constValueData(), array0.nValues());
    
    values_array += array0.nValues();
    
    quickCopy(values_array, array1.constValueData(), array1.nValues());
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const PackedArray2D<T> &other)
                 : d(other.d)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::~PackedArray2D()
{}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D<T>::operator==(const PackedArray2D<T> &other) const
{
    if (d.constData() == other.d.constData())
        return true;
        
    if (this->nArrays() != other.nArrays())
        return false;
        
    const Array *this_arrays = this->constData();
    const Array *other_arrays = other.constData();
    quint32 narrays = this->count();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        if (this_arrays[i] != other_arrays[i])
            return false;
    }
    
    return true;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D<T>::operator!=(const PackedArray2D<T> &other) const
{
    return not this->operator==(other);
}

/** Assert that the index 'i' is valid for this array

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::assertValidIndex(quint32 i) const
{
    if (i >= quint32(this->count()))
        SireBase::detail::throwPackedArray2D_invalidIndex(i, this->count());
}

/** Return a reference to the ith array in this container

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& PackedArray2D<T>::operator[](quint32 i) const
{
    this->assertValidIndex(i);
    return d->arrayData()[i];
}

/** Return a reference to the ith array in this container

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& PackedArray2D<T>::at(quint32 i) const
{
    return this->operator[](i);
}

/** Return a reference to the jth value of the ith array in 
    this container
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D<T>::operator()(quint32 i, quint32 j) const
{
    return this->at(i).at(j);
}

/** Return a modifiable reference to the jth value of the ith array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& PackedArray2D<T>::operator()(quint32 i, quint32 j)
{
    this->assertValidIndex(i);
    Array &array = d->arrayData()[i];
    
    return array[j];
}

/** Return a reference to the jth value of the ith array in 
    this container
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D<T>::at(quint32 i, quint32 j) const
{
    return this->operator()(i,j);
}

/** Return the number of arrays in this PackedArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::count() const
{
    return d->nArrays();
}

/** Return the number of arrays in this PackedArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::size() const
{
    return this->count();
}

/** Return the number of arrays in this PackedArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::nArrays() const
{
    return this->count();
}

/** Return the total number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::nValues() const
{
    return d->nValues();
}

/** Return the total number of values in the ith array

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::nValues(quint32 i) const
{
    return this->at(i).nValues();
}

/** Return whether or not this array is empty (contains no arrays) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D<T>::isEmpty() const
{
    return this->count() == 0;
}

/** Return a raw pointer to the array of arrays. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* PackedArray2D<T>::data() const
{
    return d->arrayData();
}

/** Return a raw pointer to the array of arrays. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* PackedArray2D<T>::constData() const
{
    return d->arrayData();
}

/** Detach from shared storage */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::detach()
{
    d->detach();
}

/** Return a raw pointer to the array of values in the ith array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::data(quint32 i) const
{
    return this->at(i).data();
}

/** Return a modifiable pointer to the array of values in the
    ith array. *DO NOT* delete this pointer! 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D<T>::data(quint32 i)
{
    this->assertValidIndex(i);
    return d->arrayData()[i].data();
}

/** Return a raw pointer to the array of values in the ith array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::constData(quint32 i) const
{
    return this->data(i);
}

/** Return a raw pointer to the array of all values in all arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::valueData() const
{
    return d->valueData();
}

/** Return a modifiable pointer to the array of all values in 
    all arrays. *DO NOT* delete this pointer. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D<T>::valueData()
{
    return d->valueData();
}

/** Return a raw pointer to the array of all values in all arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::constValueData() const
{
    return d->valueData();
}

/** Return the QVector that has the same contents as this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVector<T> PackedArray2D<T>::toQVector() const
{
    if (this->nValues() == 0)
        return QVector<T>();
        
    QVector<T> ret( this->nValues() );
    
    quickCopy<T>(ret.data(), this->valueData(), this->nValues());
    
    return ret;
}

/** Return the QVector that has the same contents as this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVector< QVector<T> > PackedArray2D<T>::toQVectorVector() const
{
    if (this->isEmpty())
        return QVector< QVector<T> >();
        
    QVector< QVector<T> > ret( this->count() );

    QVector<T> *ret_array = ret.data();
    const Array *arrays = this->data();
    
    for (int i=0; i<this->count(); ++i)
    {
        ret_array[i] = arrays[i].toQVector();
    }
    
    return ret;
}

/** Return a string representation of this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString PackedArray2D<T>::toString() const
{
    return Sire::toString( PackedArray2D<T>::toQVectorVector() );
}

namespace detail
{

void throwCannotConvertVariantError(const char *this_type, 
                                    const char *type_t,
                                    const QString &codeloc);

template<class T>
T default_construct()
{
    return T();
}

template<>
inline double default_construct<double>()
{
    return double(0);
}

template<>
inline qint64 default_construct<qint64>()
{
    return qint64(0);
}

template<class T>
struct VariantConverter
{
    static PackedArray2D<QVariant> output(const PackedArray2D<T> &array)
    {
        if (array.isEmpty())
            return PackedArray2D<QVariant>();
            
        int narrays = array.nArrays();
        
        QVector< QVector<QVariant> > variant(narrays);
        QVector<QVariant> *variant_array = variant.data();
        const typename PackedArray2D<T>::Array *array_array = array.constData();
        
        for (int i=0; i<narrays; ++i)
        {
            const typename PackedArray2D<T>::Array &this_array = array_array[i];
        
            int nvalues = this_array.count();
            const T *this_array_array = this_array.constData();
            
            QVector<QVariant> this_variant(nvalues);
            QVariant *this_variant_array = this_variant.data();
            
            for (int j=0; j<nvalues; ++j)
            {
                this_variant_array[j] = QVariant::fromValue<T>(this_array_array[j]);
            }
            
            variant_array[i] = this_variant;
        }
        
        return PackedArray2D<QVariant>(variant);
    }
    
    static PackedArray2D<T> input(const PackedArray2D<QVariant> &variant)
    {
        if (variant.isEmpty())
            return PackedArray2D<T>();
            
        int narrays = variant.nArrays();
        
        QVector< QVector<T> > array(narrays);
        QVector<T> *array_array = array.data();
        const typename PackedArray2D<QVariant>::Array *variant_array 
                                                            = variant.constData();
        
        for (int i=0; i<narrays; ++i)
        {
            const typename PackedArray2D<QVariant>::Array &this_variant 
                                                                = variant_array[i];
        
            int nvalues = this_variant.count();
            const QVariant *this_variant_array = this_variant.constData();
            
            QVector<T> this_array(nvalues);
            T *this_array_array = this_array.data();
            
            for (int j=0; j<nvalues; ++j)
            {
                const QVariant &value = this_variant_array[j];
             
                if (value.isNull())
                {
                    //use a default-constructed value
                    this_array_array[j] = SireBase::detail::default_construct<T>();
                }
                else if (not value.canConvert<T>())
                {
                    throwCannotConvertVariantError(value.typeName(),
                                            QMetaType::typeName( qMetaTypeId<T>() ),
                                            CODELOC );
                }
                else
                    this_array_array[j] = value.value<T>();
            }
            
            array_array[i] = this_array;
        }
        
        return PackedArray2D<T>(array);
    }
};

template<>
struct SIREBASE_EXPORT VariantConverter<QVariant>
{
    static PackedArray2D<QVariant> output(const PackedArray2D<QVariant> &array)
    {
        return array;
    }
    
    static PackedArray2D<QVariant> input(const PackedArray2D<QVariant> &variant)
    {
        return variant;
    }
};

}

/** Return this array converted into an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<QVariant> PackedArray2D<T>::toVariant() const
{
    return detail::VariantConverter<T>::output(*this);
}

/** Return a PackedArray constructed from an array of QVariants */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T> PackedArray2D<T>::fromVariant(const PackedArray2D<QVariant> &variant)
{
    return detail::VariantConverter<T>::input(variant);
}

/** Update the ith array so that it has the same contents as 'array'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::update(quint32 i, const Array &array)
{
    this->assertValidIndex(i);
    d->arrayData()[i].update(array);
}

/** Update the ith array so that it has the same contents as 'array'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::update(quint32 i, const QVector<T> &array)
{
    this->assertValidIndex(i);
    d->arrayData()[i].update(array);
}

/** Update the arrays whose indicies are in 'idxs' so that they have the
    same contents as the passed arrays
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
template<class C>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::updateAll(const C &idxs, const PackedArray2D<T> &arrays)
{
    if (idxs.count() != arrays.count())
        throw SireError::incompatible_error( QObject::tr(
                "You cannot update the arrays because the number of indicies (%1) "
                "is not equal to the number of arrays (%2). Indicies are %3.")
                    .arg(idxs.count()).arg(arrays.count())
                    .arg(Sire::toString(idxs)), CODELOC );
                
    PackedArray2D<T> other(*this);
    
    const typename PackedArray2D<T>::Array *array = arrays.constData();
    int i = 0;
    
    for (typename C::const_iterator it = idxs.constBegin();
         it != idxs.constEnd();
         ++it)
    {
        other.assertValidIndex(*it);
        other.d->arrayData()[*it].update(array[i]);
        ++i;
    }
    
    this->operator=(other);
}

/** Update the arrays whose indicies are in 'idxs' so that they have the
    same contents as the passed arrays
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::updateAll(const QVarLengthArray<int> &idxs, 
                                 const PackedArray2D<T> &arrays)
{
    if (idxs.count() != arrays.count())
        throw SireError::incompatible_error( QObject::tr(
                "You cannot update the arrays because the number of indicies (%1) "
                "is not equal to the number of arrays (%2).")
                    .arg(idxs.count()).arg(arrays.count()), CODELOC );
                
    PackedArray2D<T> other(*this);
    
    const typename PackedArray2D<T>::Array *array = arrays.constData();
    
    for (int i=0; i<idxs.count(); ++i)
    {
        int idx = idxs.constData()[i];
    
        other.assertValidIndex(idx);
        other.d->arrayData()[idx].update(array[i]);
    }
    
    this->operator=(other);
}

/** Update the arrays whose indicies are in 'idxs' so that they have the
    same contents as the passed arrays
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
template<class T>
template<class C>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::updateAll(const C &idxs, const QVector< QVector<T> > &arrays)
{
    if (idxs.count() != arrays.count())
        throw SireError::incompatible_error( QObject::tr(
                "You cannot update the arrays because the number of indicies (%1) "
                "is not equal to the number of arrays (%2). Indicies are %3.")
                    .arg(idxs.count()).arg(arrays.count())
                    .arg(Sire::toString(idxs)), CODELOC );
                
    PackedArray2D<T> other(*this);
    
    typename QVector< QVector<T> >::const_iterator array = arrays.constBegin();
    
    for (typename C::const_iterator it = idxs.constBegin();
         it != idxs.constEnd();
         ++it)
    {
        other.assertValidIndex(*it);
        other.d->arrayData()[*it].update(*array);
        ++array;
    }
    
    this->operator=(other);
}

/** Append the passed array onto the end of the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::append(const Array &array)
{
    this->operator=( PackedArray2D<T>(*this, array) );
}

/** Append the passed arrays onto the end of the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::append(const PackedArray2D<T> &arrays)
{
    this->operator=( PackedArray2D<T>(*this, arrays) );
}

/** Append the passed array onto the end of the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::append(const QVector<T> &array)
{
    this->operator=( PackedArray2D<T>(*this, array) );
}

/** Append the passed arrays onto the end of the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::append(const QVector< QVector<T> > &arrays)
{
    this->operator=( PackedArray2D<T>(*this, arrays) );
}

/** Remove the array at the specified index

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::remove(quint32 idx)
{
    quint32 nvals = this->nValues() - this->nValues(idx);
    quint32 narrays = this->nArrays() - 1;
    
    if (nvals == 0 or narrays == 0)
    {
        this->operator=( PackedArray2D<T>() );
        return;
    }
    
    detail::SharedArray2DPtr< detail::PackedArray2DData<T> > new_d
                                     = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = new_d.data();

    const typename PackedArray2D<T>::Array *array_data = this->constData();
    
    int int_idx = idx;
    
    //dimension each packed array
    for (int i=0; i<this->nArrays(); ++i)
    {
        if (i < int_idx)
            dptr->setNValuesInArray(i, array_data[i].count());
        
        else if (i > int_idx)
            dptr->setNValuesInArray(i-1, array_data[i].count());
    }
    
    dptr->close();
    
    //now copy all of the data
    T *values_array = dptr->valueData();

    for (int i=0; i<this->nArrays(); ++i)
    {
        if (i != int_idx)
        {
            quickCopy(values_array, array_data[i].constData(),
                      array_data[i].count());
            
            values_array += array_data[i].count();
        }
    }
    
    d = new_d;
}

/** Remove all of the arrays at the specified indicies

    \throw SireError::invalid_index
*/
template<class T>
template<class C>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::removeAll(const C &idxs)
{
    if (idxs.isEmpty())
        return;
    
    else if (idxs.count() == 1)
    {
        this->remove( *(idxs.constBegin()) );
        return;
    }
    
    QVector<bool> to_keep(this->nArrays(), true);
    
    for (typename C::const_iterator it = idxs.constBegin(); 
         it != idxs.constEnd();
         ++it)
    {
        this->assertValidIndex(*it);
        to_keep[*it] = false;
    }
    
    quint32 narrays = 0;
    
    for (QVector<bool>::const_iterator it = to_keep.constBegin();
         it != to_keep.constEnd();
         ++it)
    {
        if (*it)
            narrays += 1;
    }
    
    if (narrays == 0)
    {
        this->operator=( PackedArray2D<T>() );
        return;
    }
    
    quint32 nvals = 0;
    
    const typename PackedArray2D<T>::Array *array_data = this->constData();
    
    for (int i=0; i<this->nArrays(); ++i)
    {
        if (to_keep.constData()[i])
        {
            nvals += array_data[i].count();
        }
    }
    
    if (nvals == 0)
    {
        this->operator=( PackedArray2D<T>() );
        return;
    }
    
    detail::SharedArray2DPtr< detail::PackedArray2DData<T> > new_d
                                     = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = new_d.data();

    int idx = 0;
    
    //dimension each packed array
    for (int i=0; i<this->nArrays(); ++i)
    {
        if (to_keep.constData()[i])
        {
            dptr->setNValuesInArray(idx, array_data[i].count());
            ++idx;
        }
    }
    
    dptr->close();
    
    //now copy all of the data
    T *values_array = dptr->valueData();

    for (int i=0; i<this->nArrays(); ++i)
    {
        if (to_keep.constData()[i])
        {
            quickCopy(values_array, array_data[i].constData(),
                      array_data[i].count());
            
            values_array += array_data[i].count();
        }
    }
    
    d = new_d;
}

/** Remove all of the arrays at the specified indicies

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::removeAll(const QVarLengthArray<int> &idxs)
{
    if (idxs.isEmpty())
        return;
    
    else if (idxs.count() == 1)
    {
        this->remove( idxs.constData()[0] );
        return;
    }
    
    QVector<bool> to_keep(this->nArrays(), true);
    
    for (int i=0; i<idxs.count(); ++i)
    {
        this->assertValidIndex( idxs.constData()[i] );
        to_keep[ idxs.constData()[i] ] = false;
    }
    
    quint32 narrays = 0;
    
    for (QVector<bool>::const_iterator it = to_keep.constBegin();
         it != to_keep.constEnd();
         ++it)
    {
        if (*it)
            narrays += 1;
    }
    
    if (narrays == 0)
    {
        this->operator=( PackedArray2D<T>() );
        return;
    }
    
    quint32 nvals = 0;
    
    const typename PackedArray2D<T>::Array *array_data = this->constData();
    
    for (int i=0; i<this->nArrays(); ++i)
    {
        if (to_keep.constData()[i])
        {
            nvals += array_data[i].count();
        }
    }
    
    if (nvals == 0)
    {
        this->operator=( PackedArray2D<T>() );
        return;
    }
    
    detail::SharedArray2DPtr< detail::PackedArray2DData<T> > new_d
                                     = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = new_d.data();

    int idx = 0;
    
    //dimension each packed array
    for (int i=0; i<this->nArrays(); ++i)
    {
        if (to_keep.constData()[i])
        {
            dptr->setNValuesInArray(idx, array_data[i].count());
            ++idx;
        }
    }
    
    dptr->close();
    
    //now copy all of the data
    T *values_array = dptr->valueData();

    for (int i=0; i<this->nArrays(); ++i)
    {
        if (to_keep.constData()[i])
        {
            quickCopy(values_array, array_data[i].constData(),
                      array_data[i].count());
            
            values_array += array_data[i].count();
        }
    }
    
    d = new_d;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireBase

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Serialise a PackedArray2D<T>::Array to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const typename SireBase::detail::PackedArray2D_Array<T> &array)
{
    SireBase::detail::writePackedArray2DArrayHeader(ds, 1);
    
    //serialise the number of objects in this array
    ds << quint32(array.count());
    
    //serialise all of the values
    for (int i=0; i<array.count(); ++i)
    {
        ds << array.constData()[i];
    }
    
    return ds;
}

/** Extract a PackedArray2D<T>::Array from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        typename SireBase::detail::PackedArray2D_Array<T> &array)
{
    SireBase::detail::readPackedArray2DArrayHeader(ds, 1);
    
    //extract the number of objects
    quint32 nvals;
    ds >> nvals;
    
    if (nvals == 0)
    {
        array = typename SireBase::PackedArray2D<T>::Array();
        return ds;
    }
    
    //create an array of sufficient size
    typename SireBase::PackedArray2D<T>::Array new_array(nvals);
    
    //read in the values
    T *data = new_array.data();
    
    for (quint32 i=0; i<nvals; ++i)
    {
        ds >> data[i];
    }
    
    array = new_array;
    
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds,
                        const SireBase::detail::SharedArray2DPtr< 
                                    SireBase::detail::PackedArray2DData<T> > &d)
{
    SireStream::SharedDataStream sds(ds);

    quint32 narrays = d->nArrays();

    //serialise the number of arrays in this array
    sds << narrays;
    
    if (narrays == 0)
        return ds;
    
    //now go through each array and serialise the number of values
    //in each array
    for (quint32 i=0; i<narrays; ++i)
    {
        sds << quint32( d->arrayData()[i].count() );
    }
    
    //now serialise all of the data
    quint32 nvalues = d->nValues();
    const T *array_data = d->valueData();
    
    for (quint32 i=0; i<nvalues; ++i)
    {
        sds << array_data[i];
    }
    
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds,
                        SireBase::detail::SharedArray2DPtr< 
                                    SireBase::detail::PackedArray2DData<T> > &d)
{
    SireStream::SharedDataStream sds(ds);

    //get the number of arrays in this array
    quint32 narrays;
    
    sds >> narrays;
    
    if (narrays == 0)
    {
        d = SireBase::detail::getSharedNull<T>();
        return ds;
    }

    QVarLengthArray<quint32> array_sizes(narrays);

    quint32 nvals = 0;
    
    for (quint32 i=0; i<narrays; ++i)
    {
        sds >> array_sizes[i];
        nvals += array_sizes[i];
    }

    if (nvals == 0)
    {
        d = SireBase::detail::getSharedNull<T>();
        return ds;
    }
    
    SireBase::detail::SharedArray2DPtr< SireBase::detail::PackedArray2DData<T> > d2 
                                    = SireBase::detail::createArray<T>(narrays, nvals);
    
    SireBase::detail::PackedArray2DData<T> *dptr = d2.data();
    
    //dimension each packed array
    for (quint32 i=0; i<narrays; ++i)
    {
        //dimension the array
        dptr->setNValuesInArray(i, array_sizes[i]);
    }
    
    dptr->close();
    
    //now extract all of the objects...
    T *data = dptr->valueData();
    
    for (quint32 i=0; i<nvals; ++i)
    {
        sds >> data[i];
    }
    
    d = d2;
    
    return ds;
}

namespace SireBase{ namespace detail {

template<class T>
struct GetPackedArrayPointer
{
    static bool isEmpty(const SharedArray2DPtr< PackedArray2DData<T> > &d)
    {
        return d.constData() == 0;
    }

    static const void* value(const SharedArray2DPtr< PackedArray2DData<T> > &d)
    {
        return d.constData();
    }
    
    static void load(QDataStream &ds, SharedArray2DPtr< PackedArray2DData<T> > &d)
    {
        ds >> d;
    }
    
    static void save(QDataStream &ds, const SharedArray2DPtr< PackedArray2DData<T> > &d)
    {
        ds << d;
    }
};

}}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SireStream::SharedDataStream& operator>>(SireStream::SharedDataStream &sds,
                                         SireBase::detail::SharedArray2DPtr< 
                                            SireBase::detail::PackedArray2DData<T> > &d)
{
    sds.sharedLoad< SireBase::detail::SharedArray2DPtr<
                        SireBase::detail::PackedArray2DData<T> >,
                            SireBase::detail::GetPackedArrayPointer<T> >(d);

    return sds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SireStream::SharedDataStream& operator<<(SireStream::SharedDataStream &sds,
                                         const SireBase::detail::SharedArray2DPtr< 
                                            SireBase::detail::PackedArray2DData<T> > &d)
{
    sds.sharedSave< SireBase::detail::SharedArray2DPtr<
                        SireBase::detail::PackedArray2DData<T> >,
                            SireBase::detail::GetPackedArrayPointer<T> >(d);

    return sds;
}

/** Serialise a PackedArray2D<T> to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, 
                        const SireBase::PackedArray2D<T> &arrays)
{
    SireBase::detail::writePackedArray2DHeader(ds);
    
    SireStream::SharedDataStream sds(ds);
    
    sds << arrays.d;
    
    return ds;
}

/** Extract a PackedArray2D<T> from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, 
                        SireBase::PackedArray2D<T> &arrays)
{
    SireStream::VersionID v = SireBase::detail::readPackedArray2DHeader(ds);
    
    SireStream::SharedDataStream sds(ds);
    
    if (v == 2)
    {
        sds >> arrays.d;
    }
    else if (v == 1)
    {
        ds >> arrays.d;
    }
    else
        qFatal("SEVERE PROGRAM BUG IN %s", qPrintable(CODELOC));
    
    return ds;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
