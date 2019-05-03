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

#ifndef SIREVOL_COORDGROUP_H
#define SIREVOL_COORDGROUP_H

#include <QVarLengthArray>

#include "SireMaths/vector.h"

#include "SireBase/shareddatapointer.hpp"
#include "SireBase/refcountdata.h"

#include "aabox.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class CoordGroupBase;
class CoordGroup;
class CoordGroupEditor;
class CoordGroupArray;
class CoordGroupArrayArray;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::CoordGroup&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::CoordGroup&);

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::CoordGroupArray&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::CoordGroupArray&);

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::CoordGroupArrayArray&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::CoordGroupArrayArray&);

namespace SireMaths
{
class AxisSet;
class Matrix;
class Vector;
class Quaternion;
class Transform;
}

namespace SireVol
{

namespace detail
{
class CGArrayArrayData;
class CGArrayData;
class CGData;
class CGMemory;

/** This is the implicitly shared pointer class that 
    is used to hold any of the CGMemory allocated objects
    
    @author Christopher Woods
*/
template<class T>
class CGSharedPtr
{

public:
    CGSharedPtr() : ptr(0)
    {}

    CGSharedPtr(const T *p)
    {
        ptr = const_cast<T*>(p);
    
        if (ptr)
            ptr->incref();
    }
    
    CGSharedPtr(const CGSharedPtr &other) : ptr(other.ptr)
    {
        if (ptr)
            ptr->incref();
    }
    
    CGSharedPtr(CGSharedPtr &&other) : ptr(other.ptr)
    {
        other.ptr = 0;
    }
    
    ~CGSharedPtr()
    {
        if (ptr)
            ptr->decref();
    }
    
    CGSharedPtr<T>& operator=(const CGSharedPtr &other)
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
    
    CGSharedPtr<T>& operator=(CGSharedPtr &&other)
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
    
    /** Assign this pointer to point at 'weakptr' 
        but *without* changing the reference count.
        You ABSOLUTELY MUST ensure that you call 
        CGSharedPtr::weakRelease() before this 
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

};

using SireMaths::Vector;
using SireMaths::Quaternion;
using SireMaths::Matrix;
using SireMaths::AxisSet;
using SireMaths::Transform;

/** This is the base class of all CoordGroup-like classes
    (e.g. CoordGroup and CoordGroupEditor). CoordGroup classes
    hold a group of coordinates, together with an axis-aligned
    box that completely encloses all of those points. The
    class is implicitly shared, and, since it is used in the
    most performance-sensitive parts of the code, has
    a slightly more complex implicit-sharing design.

    @author Christopher Woods
*/
class SIREVOL_EXPORT CoordGroupBase
{

friend class detail::CGData; // so can see d pointer
friend class detail::CGMemory; // so can see d pointer

friend class CoordGroupArray;
friend class CoordGroupArrayArray;

public:
    ~CoordGroupBase();

    static const char *typeName()
    {
        return "SireVol::CoordGroupBase";
    }
    
    const char* what() const
    {
        return CoordGroupBase::typeName();
    }

    bool operator==(const CoordGroupBase &other) const;
    bool operator!=(const CoordGroupBase &other) const;

    QString toString() const;

    bool maybeDifferent(const CoordGroupBase &other) const;

    const Vector& at(quint32 i) const;
    const Vector& operator[](quint32 i) const;

    const AABox& aaBox() const;

    const Vector* constData() const;
    const Vector* data() const;

    bool isEmpty() const;

    int count() const;
    int size() const;

    QVector<Vector> toVector() const;

    void assertValidIndex(quint32 i) const;

    void assertSameSize(const QVector<Vector> &coordinates) const;
    void assertSameSize(const CoordGroupBase &other) const;

protected:
    CoordGroupBase();
    CoordGroupBase(detail::CGData *data);

    CoordGroupBase(quint32 size, const Vector &value = Vector());
    CoordGroupBase(quint32 size, const Vector *values);

    CoordGroupBase(const CoordGroupArray &cgarray);
    CoordGroupBase(const CoordGroupArrayArray &cgarrays);

    CoordGroupBase(const QVector<Vector> &coordinates);

    CoordGroupBase(const CoordGroupBase &other);

    CoordGroupBase& operator=(const CoordGroupBase &other);

    /** Pointer to the CGData object that describes
        this CoordGroup */
    detail::CGSharedPtr<detail::CGData> d;
};

/** This class holds a group of coordinates. This group forms the basis of the
    Molecular CutGroup, as defined in SireMol. A CoordGroup contains a list of
    coordinates, together with an AABox which provides information as to the
    center and extents of this group. SireVol is designed to calculate distances
    between points in different CoordGroups, or to calculate distances between
    points within a CoordGroup. A CoordGroup is implicitly shared and is
    designed to be fast to use, and fast to copy.

    @author Christopher Woods
*/
class SIREVOL_EXPORT CoordGroup : public CoordGroupBase
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const CoordGroup&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, CoordGroup&);

friend class CoordGroupEditor;
friend class detail::CGData; // so can see d pointer
friend class detail::CGMemory; // so can see d pointer

public:
    CoordGroup();
    CoordGroup(quint32 size);
    CoordGroup(quint32 size, const Vector &value);
    CoordGroup(quint32 size, const Vector *values);
    CoordGroup(const CoordGroupArray &cgarray);
    CoordGroup(const CoordGroupArrayArray &cgarrays);
    CoordGroup(const QVector<Vector> &points);

    CoordGroup(const CoordGroup &other);

    ~CoordGroup();

    static const char* typeName();
    
    const char* what() const
    {
        return CoordGroup::typeName();
    }

    CoordGroup& operator=(const CoordGroup &other);
    CoordGroup& operator=(CoordGroupEditor &other);

    CoordGroupEditor edit() const;

private:
    CoordGroup(const CoordGroupEditor &other);

    CoordGroup(detail::CGData *data);

    static void throwInvalidCountError(uint nats0, uint nats1);
};

/** This class is used to edit a CoordGroup. This class is used when you want to
    make several small changes to a CoordGroup, but do not want the CoordGroup to
    update its internal state after each change (e.g. you are moving each point in
    turn, and do not want the AABox to be updated for every step!)

    You use a CoordGroupEditor like this;

    \code

    //create a CoordGroup with space for 100 coordinates
    CoordGroup coordgroup(100);

    //create an editor for this group
    CoordGroupEditor editor = coordgroup.edit();

    //manipulate each coordinate in turn
    for (int i=0; i<100; ++i)
        editor[i] = Vector(i,i+1,i+2);

    //commit the changes
    coordgroup = editor.commit();

    \endcode

    @author Christopher Woods
*/
class SIREVOL_EXPORT CoordGroupEditor : public CoordGroupBase
{
public:
    CoordGroupEditor();
    CoordGroupEditor(const CoordGroup &other);

    CoordGroupEditor(const CoordGroupEditor &other);

    ~CoordGroupEditor();

    static const char* typeName();
    
    const char* what() const
    {
        return CoordGroupEditor::typeName();
    }

    CoordGroupEditor& operator=(const CoordGroup &cgroup);
    CoordGroupEditor& operator=(const CoordGroupEditor &other);

    Vector& operator[](quint32 i);

    Vector* data();

    CoordGroupEditor& translate(const Vector &delta);
    CoordGroupEditor& translate(quint32 i, const Vector &delta);

    CoordGroupEditor& rotate(const Quaternion &quat, const Vector &point);
    CoordGroupEditor& rotate(const Matrix &rotmat, const Vector &point);

    CoordGroupEditor& rotate(quint32 i, const Quaternion &quat, const Vector &point);
    CoordGroupEditor& rotate(quint32 i, const Matrix &rotmat, const Vector &point);

    CoordGroupEditor& transform(const Transform &t);
    CoordGroupEditor& transform(quint32 i, const Transform &t);

    CoordGroupEditor& setCoordinates(const QVector<Vector> &newcoords);
    CoordGroupEditor& setCoordinates(const CoordGroupBase &newcoords);

    CoordGroupEditor& setCoordinates(quint32 i, const Vector &newcoords);

    CoordGroupEditor& mapInto(const SireMaths::AxisSet &axes);
    CoordGroupEditor& mapInto(quint32 i, const SireMaths::AxisSet &axes);
    
    CoordGroupEditor& changeFrame(const SireMaths::AxisSet &from_frame,
                                  const SireMaths::AxisSet &to_frame);
                                  
    CoordGroupEditor& changeFrame(quint32 i,
                                  const SireMaths::AxisSet &from_frame,
                                  const SireMaths::AxisSet &to_frame);

    CoordGroup commit();

    operator CoordGroup()
    {
        return this->commit();
    }

private:
    /** Whether or not the AABox needs to be recalculated */
    bool needsupdate;
};

/** This class holds an array of CoordGroups. While you could  
    of course just use a QVector<CoordGroup>, this array
    optimises the memory layout of all of the CoordGroups
    so that they all lie contiguously along the same piece
    of memory (and indeed, all of the AABoxes are grouped
    together, while all of the coordinates are grouped together).
    
    The memory packing means that this array is much more
    limited than a QVector<CoordGroup>, i.e. you can't
    add or remove CoordGroups from the array, and you can't
    do anything to the contained CoordGroups except for 
    change their coordinates.
  
    This class is really meant to be used as a fast container
    that allow rapid iteration over all of the contained
    CoordGroups / coordinates
        
    @author Christopher Woods
*/
class SIREVOL_EXPORT CoordGroupArray
{

friend class detail::CGArrayData; // so can see d pointer
friend class detail::CGMemory; // so can see d pointer

friend class CoordGroupArrayArray;

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const CoordGroupArray&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, CoordGroupArray&);

public:
    CoordGroupArray();
    CoordGroupArray(const CoordGroup &cgroup);
    
    CoordGroupArray(const QVector< QVector<Vector> > &points);
    CoordGroupArray(const QVector<CoordGroup> &cgroups);
    
    CoordGroupArray(const CoordGroupArray &array0,
                    const CoordGroupArray &array1);
    
    CoordGroupArray(const CoordGroupArray &other);
    
    ~CoordGroupArray();
    
    CoordGroupArray& operator=(const CoordGroupArray &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return CoordGroupArray::typeName();
    }
    
    bool operator==(const CoordGroupArray &other) const;
    bool operator!=(const CoordGroupArray &other) const;

    const CoordGroup& operator[](quint32 i) const;
    const CoordGroup& at(quint32 i) const;

    QString toString() const;

    int count() const;
    int size() const;
    
    bool isEmpty() const;

    int nCoordGroups() const;
    int nCoords() const;

    AABox aaBox() const;

    const CoordGroup* data() const;
    const CoordGroup* constData() const;

    const Vector* coordsData() const;
    const Vector* constCoordsData() const;

    const AABox* aaBoxData() const;
    const AABox* constAABoxData() const;

    CoordGroup merge() const;
    
    void append(const CoordGroup &cgroup);
    void append(const CoordGroupArray &cgroups);
    
    void remove(quint32 i);
    void remove(quint32 i, int count);
    void remove(const QVarLengthArray<quint32> &idxs);
    
    void update(quint32 i, const CoordGroup &cgroup);
    void update(quint32 i, const QVector<Vector> &coords);
    void update(quint32 i, const Vector *coords, int ncoords);

    void translate(const Vector &delta);
    void translate(quint32 i, const Vector &delta);
    
    void rotate(const Quaternion &quat, const Vector &point);
    void rotate(const Matrix &rotmat, const Vector &point);
    
    void rotate(quint32 i, const Quaternion &quat, const Vector &point);
    void rotate(quint32 i, const Matrix &rotmat, const Vector &point);
    
    void transform(const Transform &t);
    void transform(quint32 i, const Transform &t);
    
    void mapInto(const AxisSet &axes);
    void mapInto(quint32 i, const AxisSet &axes);
    
    void changeFrame(const AxisSet &from_frame, const AxisSet &to_frame);
    void changeFrame(quint32 i,
                     const AxisSet &from_frame, const AxisSet &to_frame);

    void assertValidIndex(quint32 i) const;
    
    void assertValidCoordGroup(quint32 i) const;
    void assertValidCoordinate(quint32 i) const;

protected:
    CoordGroupArray(detail::CGArrayData *data);

    /** Implicitly shared pointer to the data for this array */
    detail::CGSharedPtr<detail::CGArrayData> d;

private:
    void pvt_remove(const QVarLengthArray<bool> &to_remove);
};

/** This class holds an array of CoordGroupArrays. This is 
    used to pack all of the CoordGroupArrays (and thus
    all contained CoordGroups) into a single contiguous
    block of memory. This should improve efficiency of 
    iterating over these groups/coordinates, but it does
    make this array less flexible than a simple
    QVector<CoordGroupArray>, or QVector< QVector<CoordGroup> >.
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT CoordGroupArrayArray
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const CoordGroupArrayArray&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, CoordGroupArrayArray&);

public:
    CoordGroupArrayArray();
    
    CoordGroupArrayArray(const CoordGroup &cgroup);
    CoordGroupArrayArray(const CoordGroupArray &cgarray);
    
    CoordGroupArrayArray(const QVector<CoordGroupArray> &cgarrays);
    CoordGroupArrayArray(const QVector< QVector<CoordGroup> > &cgarrays);
    CoordGroupArrayArray(const QVector< QVector< QVector<Vector> > > &points);
    
    CoordGroupArrayArray(const CoordGroupArrayArray &other);
    
    ~CoordGroupArrayArray();
    
    static const char* typeName();
    
    const char* what() const
    {
        return CoordGroupArrayArray::typeName();
    }
    
    CoordGroupArrayArray& operator=(const CoordGroupArrayArray &other);
    
    bool operator==(const CoordGroupArrayArray &other) const;
    bool operator!=(const CoordGroupArrayArray &other) const;

    QString toString() const;

    const CoordGroupArray& operator[](quint32 i) const;
    const CoordGroupArray& at(quint32 i) const;

    int count() const;
    int size() const;

    int nCoordGroupArrays() const;
    int nCoordGroups() const;
    int nCoords() const;

    CoordGroup merge() const;

    AABox aaBox() const;

    const CoordGroupArray* data() const;
    const CoordGroupArray* constData() const;

    const CoordGroup* coordGroupData() const;
    const CoordGroup* constCoordGroupData() const;

    const Vector* coordsData() const;
    const Vector* constCoordsData() const;

    const AABox* aaBoxData() const;
    const AABox* constAABoxData() const;
    
    void update(quint32 i, const CoordGroupArray &array);
    void update(quint32 i, quint32 j, const CoordGroup &cgroup);

    void translate(const Vector &delta);
    void translate(quint32 i, const Vector &delta);
    void translate(quint32 i, quint32 j, const Vector &delta);
    
    void rotate(const Quaternion &quat, const Vector &point);
    void rotate(const Matrix &rotmat, const Vector &point);
    
    void rotate(quint32 i, const Quaternion &quat, const Vector &point);
    void rotate(quint32 i, const Matrix &rotmat, const Vector &point);

    void rotate(quint32 i, quint32 j,
                const Quaternion &quat, const Vector &point);
    void rotate(quint32 i, quint32 j,
                const Matrix &rotmat, const Vector &point);
    
    void transform(const Transform &t);
    void transform(quint32 i, const Transform &t);
    void transform(quint32 i, quint32 j, const Transform &t);
    
    void mapInto(const AxisSet &axes);
    void mapInto(quint32 i, const AxisSet &axes);
    void mapInto(quint32 i, quint32 j, const AxisSet &axes);
    
    void changeFrame(const AxisSet &from_frame, const AxisSet &to_frame);
    void changeFrame(quint32 i,
                     const AxisSet &from_frame, const AxisSet &to_frame);
    void changeFrame(quint32 i, quint32 j,
                     const AxisSet &from_frame, const AxisSet &to_frame);

    void assertValidIndex(quint32 i) const;
    
    void assertValidCoordGroupArray(quint32 i) const;
    void assertValidCoordGroup(quint32 i) const;
    void assertValidCoordinate(quint32 i) const;
    
    void assertValidCoordGroup(quint32 i, quint32 j) const;
    
private:
    /** Implicitly shared pointer to the array data */
    detail::CGSharedPtr<detail::CGArrayArrayData> d;
};

}

Q_DECLARE_METATYPE( SireVol::CoordGroup );
Q_DECLARE_METATYPE( SireVol::CoordGroupEditor );
Q_DECLARE_METATYPE( SireVol::CoordGroupArray );
Q_DECLARE_METATYPE( SireVol::CoordGroupArrayArray );

SIRE_EXPOSE_CLASS( SireVol::CoordGroupBase )
SIRE_EXPOSE_CLASS( SireVol::CoordGroup )
SIRE_EXPOSE_CLASS( SireVol::CoordGroupEditor )
SIRE_EXPOSE_CLASS( SireVol::CoordGroupArray )
SIRE_EXPOSE_CLASS( SireVol::CoordGroupArrayArray )

SIRE_END_HEADER

#endif
