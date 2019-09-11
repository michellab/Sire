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

#include "coordgroup.h"

#include "SireMaths/matrix.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/axisset.h"
#include "SireMaths/rotate.h"
#include "SireMaths/align.h"

#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireVol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

namespace SireVol
{
namespace detail
{

class CGMemory;

class CGData;
class CGArrayData;
class CGArrayArrayData;

/** This class provides the main metadata that describes the array
    of array of CoordGroups */
class CGArrayArrayData : public RefCountData
{
friend class CGMemory;

public:
    CGArrayArrayData();
    CGArrayArrayData(quint32 narray, quint32 ncgroups, quint32 ncoords);
    
    CGArrayArrayData(const CGArrayArrayData &other);
    
    ~CGArrayArrayData();
    
    void incref();
    void decref();
    
    CGArrayArrayData* detach();
    
    const CGArrayData* nullArray() const;
    const CGData* nullCGroup() const;
    
    const CGData* cgDataData() const;
    const CGArrayData* cgArrayDataData() const;
    
    const CoordGroupArray* cgArrayData() const;
    const CoordGroup* cGroupData() const;
    const AABox* aaBoxData() const;
    const Vector* coordsData() const;

    CGData* cgDataData();
    CGArrayData* cgArrayDataData();
    
    CoordGroupArray* cgArrayData();
    CoordGroup* cGroupData();
    AABox* aaBoxData();
    Vector* coordsData();
    
    void setNCGroupsInArray(quint32 i, quint32 ncgroups);
    void setNPointsInCGroup(quint32 i, quint32 npoints);
    
    void close();
    
    quint32 nCGArrays() const;
    quint32 nCGroups() const;
    quint32 nCoords() const;
    
    char *memory();
    const char* memory() const;
    
private:
    /** The index in the storage array of the first CoordGroupArray  
        in this array */
    quint32 cgarray0;
    
    /** The index in the storage array of the first CGArrayData */
    quint32 cgarraydata0;
    
    /** The number of CoordGroupArrays in this array */
    quint32 ncgarrays;
    
    /** The index in the storage array of the first CoordGroup in 
        this array */
    quint32 cgroup0;
    
    /** The index in the storage array for the first CGData in this array */
    quint32 cgdata0;
    
    /** The index in the storage array of the first AABox in this
        array (there is one for each CoordGroup) */
    quint32 aabox0;
    
    /** The number of CoordGroups in this array */
    quint32 ncgroups;
    
    /** The index in the storage array of the first vector in
        this array */
    quint32 coords0;
    
    /** The number of vectors in this array */
    quint32 ncoords;
};

/** This class provides the metainformation that describes an
    array of CoordGroups */
class CGArrayData
{

friend class CGMemory;
friend class CGArrayArrayData;

public:
    CGArrayData();
    CGArrayData(quint32 this_idx);
    
    CGArrayData(const CGArrayData &other);
    
    ~CGArrayData();

    void incref();
    void decref();
    
    CGArrayData* detach();
    
    CGArrayArrayData* extract() const;
    
    const CoordGroup* cGroupData() const;
    const AABox* aaBoxData() const;
    const Vector* coordsData() const;
    
    const CoordGroup* constCGroupData() const;
    const AABox* constAABoxData() const;
    const Vector* constCoordsData() const;
    
    CoordGroup* cGroupData();
    AABox* aaBoxData();
    Vector* coordsData();
    
    quint32 nCGroups() const;
    quint32 nCoords() const;
    
    char* memory();
    const char* memory() const;
    
private:
    /** The index in the storage array of this CoordGroupArray */
    quint32 this_cgarray;
    
    /** The index into the storage array of the first CoordGroup
        of this array */
    quint32 cgroup0;
    
    /** The index into the storage array of the first AABox of
        this array (same number as number of CoordGroups) */
    quint32 aabox0;
    
    /** The number of CoordGroups in this array */
    quint32 ncgroups;
    
    /** The index into the storage array of the first vector in 
        this array */
    quint32 coords0;
    
    /** The number of vectors in this array */
    quint32 ncoords;
};

/** This class provides the metainformation that describes
    a single CoordGroup */
class CGData
{

friend class CGMemory;
friend class CGArrayArrayData;

public:
    CGData();
    CGData(quint32 this_idx);
    
    CGData(const CGData &other);
    
    ~CGData();
    
    void incref();
    void decref();
    
    CGData* detach();
    
    CGArrayArrayData* extract() const;
    
    const char* memory() const;
    const Vector* coordsData() const;
    const AABox* aaBox() const;

    char* memory();
    Vector* coordsData();
    AABox* aaBox();
    
    quint32 nCoords() const;
    
private:
    /** The index in the storage array of this CoordGroup */
    quint32 this_cgroup;
    
    /** The index in the storage array of the first vector in this CoordGroup */
    quint32 coords0;
    
    /** The number of vectors in this CoordGroup */
    quint32 ncoords;
    
    /** The index in the storage array of the AABox for this CoordGroup */
    quint32 aabox;
};

/** This class is used to create, copy and manage the memory of 
    the CoordGroup-type classes
    
    @author Christopher Woods
*/
class CGMemory
{
public:
    static char* create(quint32 narray, quint32 ncgroups, quint32 ncoords);

    static char* detach(char *this_ptr, quint32 this_idx);

    static void destroy(CGArrayArrayData *array);
    
    static void incref(char *this_ptr, quint32 this_idx);
    static void decref(char *this_ptr, quint32 this_idx);

    static char* getRoot(char *this_ptr, quint32 this_idx);
    static const char* getRoot(const char *this_ptr, quint32 this_idx);

private:
    static quint32 getSize(quint32 narrays, quint32 ngroups, quint32 ncoords);
};

} // end of namespace detail;
} // end of namespace SireVol;

using namespace SireVol::detail;

////////
//////// Implementation of CGMemory
////////

/** Private function used to calculate the size of the memory required
    for a container with narrays CoordGroupArrays, containing in total
    ngroups CoordGroups and ncoords Vector coordinates */
quint32 CGMemory::getSize(quint32 narrays, quint32 ngroups, quint32 ncoords)
{
    //need to add space for null arrays and groups
    if (narrays == 0)
        narrays = 1;
        
    if (ngroups == 0)
        ngroups = 1;

    return sizeof(CGArrayArrayData) +
                  narrays * sizeof(CoordGroupArray) + 
                  narrays * sizeof(CGArrayData) + 
                  ngroups * sizeof(CoordGroup) +
                  ngroups * sizeof(CGData) + 
                  ngroups * sizeof(AABox) + 
                  ncoords * sizeof(Vector) + 
                  16; // allow 16 extra bytes for spacing
}

/** Create space for a container that holds narrays CoordGroupArrays,
    ncgroups CoordGroups (in total) and ncoords Vectors (in total) 
    
    // The whole point of this set of classes is so that all of the 
    // relavant coordinate data of a set of molecules is all packed
    // together in memory. This should help keep the cache full,
    // as we won't have to keep bouncing randomly around main memory
    // to find the coordinates of the next CoordGroup to test.
    //
    // The layout in memory for this data is as follows;
    //
    // ---------------------------------------------------------------------
    // |A|  B  |    C      |    D     |                E                   |
    // ---------------------------------------------------------------------
    //
    //  A = CGArrayArrayData object       ( 1 * sizeof(CGArrayArrayData) )
    //  B = Array of CoordGroupArrays +   ( narrays * sizeof(CoordGroupArray) )
    //      Array of CGArrayData objects  ( narrays * sizeof(CGArrayData) )
    //  C = Array of CoordGroups +        ( ncgroups * sizeof(CoordGroup) )
    //      Array of CGData objects       ( ncgroups * sizeof(CGData) )
    //  D = Array of AABoxes 
    //       (one for each CoordGroup)    ( ncgroups * sizeof(AABox) )
    //  E = Array of Vectors              ( ncoords * sizeof(Vector) )
    //
    // This layout packs all of the coordinates together, all of the
    // AABoxes together (used for sphere cutting tests), and then
    // all of the CoordGroup, then CGArray metadata together. This layout
    // also then allows easy iteration over sub-parts, e.g. iterating
    // over all coordinates in all CoordGroups in all arrays, or
    // over all coordinates in all CoordGroups in a specific array,
    // or over all coordinates in a specific CoordGroup.
*/
char* CGMemory::create(quint32 narrays, quint32 ncgroups, quint32 ncoords)
{
    //calculate how big the memory requirements are
    quint32 sz = getSize(narrays, ncgroups, ncoords);
    
    //allocate this space
    char *storage = new char[sz];
    
    quint32 spare_space = 16;
    
    try
    {
        //now we need to construct each object in turn at its correct
        //location in this space
        CGArrayArrayData *cgarrayarray 
                    = new (storage) CGArrayArrayData(narrays, ncgroups, ncoords);
                    
        //advance the index into the storage array so that we are now
        //pointing just after the CGArrayArrayData
        quint32 idx = sizeof(CGArrayArrayData);
                    
        //this is the location of the first CoordGroupArray :-)
        cgarrayarray->cgarray0 = idx;
        
        //the first CGArrayData lies at this index + narrays*sizeof(CoordGroupArray)
        if (narrays > 0)
            cgarrayarray->cgarraydata0 = idx + (narrays * sizeof(CoordGroupArray));

        if (narrays > 0)
        {
            quint32 dataidx = cgarrayarray->cgarraydata0;
        
            //loop over each CGArray and create it in its place
            for (quint32 i=0; i<narrays; ++i)
            {   
                //assert that there is sufficient space in the array
                BOOST_ASSERT(idx + sizeof(CoordGroupArray) < sz);
                BOOST_ASSERT(dataidx + sizeof(CGArrayData) < sz);
            
                //create the CGArrayData, letting it know where it
                //is relative to the beginning of the storage array
                CGArrayData *data = new (storage + dataidx) CGArrayData(dataidx);
                
                //now create the CoordGroupArray that uses this data
                new (storage + idx) CoordGroupArray(data);
        
                //advance the index into the storage array to point
                //just after the just-created CGArrayData
                idx += sizeof(CoordGroupArray);
                dataidx += sizeof(CGArrayData);
            }
            
            //make idx point just after all of the CGArrayDatas...
            BOOST_ASSERT( idx == cgarrayarray->cgarray0 
                                + narrays*sizeof(CoordGroupArray) );
                                
            idx = dataidx;
        }
        else
        {
            //we need to create space for the null CoordGroupArray
            quint32 data_idx = idx + sizeof(CoordGroupArray);
            CGArrayData *cgarray = new (storage + data_idx) CGArrayData(data_idx);
            
            //tell the parent where this null data is located
            cgarrayarray->cgarraydata0 = data_idx;
            
            new (storage + idx) CoordGroupArray(cgarray);
            
            cgarray->ncgroups = 0;
            cgarray->cgroup0 = 0;
            cgarray->ncoords = 0;
            cgarray->coords0 = 0;
            cgarray->aabox0 = 0;

            cgarrayarray->ncgarrays = 0;
            
            idx += sizeof(CGArrayData) + sizeof(CoordGroupArray);
        }
        
        //we are now at the location of the first CGData
        cgarrayarray->cgroup0 = idx;
        cgarrayarray->cgdata0 = idx + (ncgroups * sizeof(CoordGroup));
        
        if (ncgroups > 0)
        {
            quint32 dataidx = cgarrayarray->cgdata0;
        
            //loop over each CoordGroup and create it in its place
            for (quint32 i=0; i<ncgroups; ++i)
            {
                //assert that there is sufficient space in the array
                BOOST_ASSERT(idx + sizeof(CoordGroup) < sz);
                BOOST_ASSERT(dataidx + sizeof(CGData) < sz);
            
                //create the CGData, letting it know where it is
                //relative to the beginning of the storage array
                CGData *data = new (storage + dataidx) CGData(dataidx);
                
                //create the CoordGroup
                new (storage + idx) CoordGroup(data);
                
                //advance the index into the storage array to point
                //just after the just-created CGData
                idx += sizeof(CoordGroup);
                dataidx += sizeof(CGData); 
            }
        
            BOOST_ASSERT( idx == cgarrayarray->cgroup0 + 
                                    ncgroups * sizeof(CoordGroup) );
        
            idx = dataidx;
        
            //we are now at the location of the first AABox
            cgarrayarray->aabox0 = idx;
        
            //loop over each AABox and create it in place
            for (quint32 i=0; i<ncgroups; ++i)
            {
                BOOST_ASSERT(idx + sizeof(AABox) < sz);
            
                new (storage + idx) AABox();
            
                idx += sizeof(AABox);
            }
        
            //we are now at the location of the first Vector
            idx += 16; // seems to like it best if there is 16 bytes of space
            spare_space -= 16;
            
            cgarrayarray->coords0 = idx;
        
            //loop over each vector and create it in place
            for (quint32 i=0; i<ncoords; ++i)
            {
                BOOST_ASSERT(idx + sizeof(Vector) <= sz);

                new (storage + idx) Vector(0.0);
                
                idx += sizeof(Vector);
            }
        }
        else
        {
            //create space for the null CoordGroup
            quint32 dataidx = idx + sizeof(CoordGroup);
            
            CGData *cgroup = new (storage + dataidx) CGData(dataidx);
            new (storage + idx) CoordGroup(cgroup);
            
            //tell the parent where this null CoordGroup is located
            cgarrayarray->cgdata0 = dataidx;
            
            cgroup->ncoords = 0;
            cgroup->coords0 = 0;

            idx += sizeof(CGData) + sizeof(CoordGroup);
            
            cgroup->aabox = idx;
            
            new (storage + idx) AABox();
            
            cgarrayarray->aabox0 = idx;
            cgarrayarray->coords0 = 0;
            cgarrayarray->ncgroups = 0;
            
            idx += sizeof(AABox);
        }
        
        //we should now be at the end of the storage
        BOOST_ASSERT( idx + spare_space == sz );
                                      
        return storage;
    } 
    catch(...)
    {
        //release the memory so we don't leak
        delete[] storage;
        throw;
    }
    
    return 0;
}

/** This converts the pointer to the object
    that starts at 'this_ptr' to a pointer to the first element
    of the storage array (given the location of the object at
    index 'this_idx' in the storage array) */
char* CGMemory::getRoot(char *this_ptr, quint32 this_idx)
{
    return this_ptr - this_idx;
}

/** This converts the pointer to the object
    that starts at 'this_ptr' to a pointer to the first element
    of the storage array (given the location of the object at
    index 'this_idx' in the storage array) */
const char* CGMemory::getRoot(const char *this_ptr, quint32 this_idx)
{
    return this_ptr - this_idx;
}

/** Destroy the object 'array' */
void CGMemory::destroy(CGArrayArrayData *array)
{
    //we need to delete it in the opposite order to creation - first
    //lets delete all of the vectors
    quint32 ncoords = array->nCoords();

    char *storage = (char*)array;
        
    if (ncoords > 0)
    {
        Vector *coords = array->coordsData();
        
        for (qint32 i=ncoords-1; i>=0; --i)
        {
            coords[i].~Vector();
        }
    }
        
    //now delete all of the AABoxes
    quint32 ncgroups = array->nCGroups();
        
    if (ncgroups > 0)
    {
        AABox *aaboxes = array->aaBoxData();
            
        for (qint32 i=ncgroups-1; i>=0; --i)
        {
            aaboxes[i].~AABox();
        }
        
        //now delete all of the CGData objects
        CoordGroup *cgroups = array->cGroupData();
            
        for (qint32 i=ncgroups-1; i>=0; --i)
        {
            CGData *data = const_cast<CGData*>( cgroups[i].d.constData() );
            
            //delete the CGData
            data->~CGData();
            
            //remove the shared pointer
            cgroups[i].d.weakRelease();
            
            //delete the CoordGroup
            cgroups[i].~CoordGroup();
        }
    }
        
    //now delete all of the CGArrayData objects
    quint32 narrays = array->nCGArrays();
        
    if (narrays > 0)
    {
        CoordGroupArray *arrays = array->cgArrayData();
            
        for (qint32 i=narrays-1; i>=0; --i)
        {
            CGArrayData *data = const_cast<CGArrayData*>(arrays[i].d.constData());
            
            //delete the data
            data->~CGArrayData();
            
            //release the shared pointer
            arrays[i].d.weakRelease();
            
            //delete the array
            arrays[i].~CoordGroupArray();
        }
    }
        
    //delete the CGArrayArrayData object itself
    array->~CGArrayArrayData();
        
    //finally delete the memory itself
    delete[] storage;
}

/** Detach this object from shared storage and return a new pointer to it. */
char* CGMemory::detach(char *this_ptr, quint32 this_idx)
{
    //get a pointer to the start of the storage for this container
    char *storage = getRoot(this_ptr, this_idx);
    
    //The CGArrayArrayData object is at the beginning of this storage array
    CGArrayArrayData *cgarrayarray = (CGArrayArrayData*) storage;

    if (cgarrayarray->ref.hasSingleReference())
    {
        //only one reference, so no need to clone
        return this_ptr;
    }
    else
    {
        //there is more than one reference to this data - it will have to 
        //be cloned - get the size of memory to be cloned
        int sz = getSize(cgarrayarray->nCGArrays(), 
                         cgarrayarray->nCGroups(),
                         cgarrayarray->nCoords());
                         
        //create space for the copy
        char *new_storage = new char[sz];
        
        //copy the data
        char *output = quickCopy<char>( new_storage, storage, sz );
        BOOST_ASSERT(output == new_storage);

        //the first part of the data is the CGArrayArrayData object
        CGArrayArrayData *new_cgarrayarray = (CGArrayArrayData*) new_storage;

        //set the reference count of this copy to 1 - note that the quickCopy
        //has copied the reference count, so this has to be reset to 1
        new_cgarrayarray->ref.reset();
        new_cgarrayarray->ref.ref();

        //now loose a reference to the original
        CGMemory::decref(this_ptr, this_idx);
        
        //the final step is to update all of the CoordGroup and CoordGroupArray
        //pointers that exist in this array (otherwise they will all
        //point to the original!)
        
        if (new_cgarrayarray->nCGArrays() > 0)
        {
            CoordGroupArray *cgarray = new_cgarrayarray->cgArrayData();
            CGArrayData *data = new_cgarrayarray->cgArrayDataData();
        
            for (quint32 i=0; i < new_cgarrayarray->nCGArrays(); ++i)
            {
                cgarray[i].d.weakRelease();
                cgarray[i].d.weakAssign( &(data[i]) );
            }
        }
        else
        {
            //update the null CoordGroupArray
            CoordGroupArray *cgarray = new_cgarrayarray->cgArrayData();
            CGArrayData *data = new_cgarrayarray->cgArrayDataData();

            cgarray[0].d.weakRelease();
            cgarray[0].d.weakAssign( &(data[0]) );
        }
        
        if (new_cgarrayarray->nCGroups() > 0)
        {
            CoordGroup *cgroups = new_cgarrayarray->cGroupData();
            CGData *data = new_cgarrayarray->cgDataData();
            
            for (quint32 i=0; i < new_cgarrayarray->nCGroups(); ++i)
            {
                cgroups[i].d.weakRelease();
                cgroups[i].d.weakAssign( &(data[i]) );
            }
        }
        else
        {
            //update the null CoordGroup
            CoordGroup *cgroups = new_cgarrayarray->cGroupData();
            CGData *data = new_cgarrayarray->cgDataData();

            cgroups[0].d.weakRelease();
            cgroups[0].d.weakAssign( &(data[0]) );
        }

        //return a pointer to the clone
        return new_storage + this_idx;
    }
}

/** Increase the reference count */
void CGMemory::incref(char *this_ptr, quint32 this_idx)
{
    ( (CGArrayArrayData*)(getRoot(this_ptr,this_idx)) )->incref();
}

/** Decrease the reference count */
void CGMemory::decref(char *this_ptr, quint32 this_idx)
{
    ( (CGArrayArrayData*)(getRoot(this_ptr,this_idx)) )->decref();
}

////////
//////// Implementation of CGArrayArrayData
////////

/** Null constructor */
CGArrayArrayData::CGArrayArrayData() 
                 : RefCountData(),
                   cgarray0(0), cgarraydata0(0), ncgarrays(0),
                   cgroup0(0), cgdata0(0), aabox0(0), ncgroups(0),
                   coords0(0), ncoords(0)
{}

/** Construct metadata for the specified number of subgroups */
CGArrayArrayData::CGArrayArrayData(quint32 narray, quint32 ngroup, 
                                   quint32 ncoord)
                 : RefCountData(),
                   cgarray0(0), cgarraydata0(0), ncgarrays(narray),
                   cgroup0(0), cgdata0(0), aabox0(0), ncgroups(ngroup),
                   coords0(0), ncoords(ncoord)
{}
  
/** Copy constructor */               
CGArrayArrayData::CGArrayArrayData(const CGArrayArrayData &other)
                 : RefCountData(),
                   cgarray0(other.cgarray0), cgarraydata0(other.cgarraydata0),
                   ncgarrays(other.ncgarrays),
                   cgroup0(other.cgroup0), cgdata0(other.cgdata0),
                   aabox0(other.aabox0),
                   ncgroups(other.ncgroups),
                   coords0(other.coords0), ncoords(other.ncoords)
{}

/** Destructor */
CGArrayArrayData::~CGArrayArrayData()
{}

/** Return a pointer to the start of the memory storage for this container */
const char* CGArrayArrayData::memory() const
{
    return CGMemory::getRoot( (const char*)this, 0 );
}

/** Return a modifiable pointer to the start of the memory storage
    for this container - this detaches this container from implicitly
    shared storage */
char* CGArrayArrayData::memory()
{
    return CGMemory::getRoot( (char*)this, 0 );
}

/** Return a pointer to the null CoordGroupArray */
const CGArrayData* CGArrayArrayData::nullArray() const
{
    BOOST_ASSERT( ncgarrays == 0 );
    BOOST_ASSERT( cgarray0 != 0 );
    
    return (const CGArrayData*)( &(memory()[cgarraydata0]) );
}

/** Return a pointer to the null CoordGroup */
const CGData* CGArrayArrayData::nullCGroup() const
{
    BOOST_ASSERT( ncgroups == 0 );
    BOOST_ASSERT( cgdata0 != 0 );
    
    return (const CGData*)( &(memory()[cgdata0]) );
}

/** Return a pointer to a CGArrayArrayData copy of this that has
    been detached from shared storage */
CGArrayArrayData* CGArrayArrayData::detach()
{
    return (CGArrayArrayData*)( CGMemory::detach( (char*)(this), 0 ) );
}

/** Increment the reference count for this object */
void CGArrayArrayData::incref()
{
    this->ref.ref();
}
    
/** Decrease the reference count for this object - this will
    delete this object if the reference count drops to zero! */
void CGArrayArrayData::decref()
{
    if (not this->ref.deref())
    {
        CGMemory::destroy(this);
    }
}

/** Return a pointer to the first CGArrayData in this container. This
    returns 0 if there are no arrays in this container */
CGArrayData* CGArrayArrayData::cgArrayDataData()
{
    if (cgarraydata0 == 0)
        return 0;
    else
        return (CGArrayData*)( memory() + cgarraydata0 );
}

/** Return a pointer to the first CGArrayData in this container. This
    returns 0 if there are no arrays in this container */
const CGArrayData* CGArrayArrayData::cgArrayDataData() const
{
    if (cgarraydata0 == 0)
        return 0;
    else
        return (const CGArrayData*)( memory() + cgarraydata0 );
}

/** Return a pointer to the first CoordGroupArray in this container. This returns
    0 if there are no arrays in this container */
const CoordGroupArray* CGArrayArrayData::cgArrayData() const
{
    if (cgarray0 == 0)
        return 0;
    else
        return (const CoordGroupArray*)( memory() + cgarray0 );
}

/** Return a pointer to the first CGData in this container. This returns 
    0 if there are no CoordGroups in this container */
CGData* CGArrayArrayData::cgDataData()
{
    if (cgdata0 == 0)
        return 0;
    else
        return (CGData*)( memory() + cgdata0 );
}

/** Return a pointer to the first CGData in this container. This returns 
    0 if there are no CoordGroups in this container */
const CGData* CGArrayArrayData::cgDataData() const
{
    if (cgdata0 == 0)
        return 0;
    else
        return (const CGData*)( memory() + cgdata0 );
}

/** Return a pointer to the first CoordGroup in this container. This returns
    0 if there are no CoordGroups in this container */
const CoordGroup* CGArrayArrayData::cGroupData() const
{
    if (cgroup0 == 0)
        return 0;
    else
        return (const CoordGroup*)( memory() + cgroup0 );
}

/** Return a pointer to the first AABox in this container. This returns
    0 if there are no AABoxes in this container */
const AABox* CGArrayArrayData::aaBoxData() const
{
    if (aabox0 == 0)
        return 0;
    else
        return (const AABox*)( memory() + aabox0 );
}

/** Return a pointer to the first Vector in this container. This
    returns 0 if there are no vectors in this container */
const Vector* CGArrayArrayData::coordsData() const
{
    if (coords0 == 0)
        return 0;
    else
        return (const Vector*)( memory() + coords0 );
}

/** Return a pointer to the first CoordGroupArray in this container. This returns
    0 if there are no arrays in this container */
CoordGroupArray* CGArrayArrayData::cgArrayData()
{
    if (cgarray0 == 0)
        return 0;
    else
        return (CoordGroupArray*)( memory() + cgarray0 );
}

/** Return a pointer to the first CoordGroup in this container. This returns
    0 if there are no CoordGroups in this container */
CoordGroup* CGArrayArrayData::cGroupData()
{
    if (cgroup0 == 0)
        return 0;
    else
        return (CoordGroup*)( memory() + cgroup0 );
}

/** Return a pointer to the first AABox in this container. This returns
    0 if there are no AABoxes in this container */
AABox* CGArrayArrayData::aaBoxData()
{
    if (aabox0 == 0)
        return 0;
    else
        return (AABox*)( memory() + aabox0 );
}

/** Return a pointer to the first Vector in this container. This
    returns 0 if there are no vectors in this container */
Vector* CGArrayArrayData::coordsData()
{
    if (coords0 == 0)
        return 0;
    else
        return (Vector*)( memory() + coords0 );
}

/** Return the number of CoordGroup arrays in this container */
quint32 CGArrayArrayData::nCGArrays() const
{
    return ncgarrays;
}

/** Return the number of CoordGroups in this container (in total) */
quint32 CGArrayArrayData::nCGroups() const
{
    return ncgroups;
}

/** Return the total number of coordinates in this container */
quint32 CGArrayArrayData::nCoords() const
{
    return ncoords;
}

/** Set the number of CoordGroups in the ith CoordGroupArray */
void CGArrayArrayData::setNCGroupsInArray(quint32 i, quint32 ncgroups)
{
    BOOST_ASSERT( i < ncgarrays );
    BOOST_ASSERT( this->cgArrayDataData() != 0 );
    
    this->cgArrayDataData()[i].ncgroups = ncgroups;
}

/** Set the number of coordinates in the ith CoordGroup */
void CGArrayArrayData::setNPointsInCGroup(quint32 i, quint32 ncoords)
{
    BOOST_ASSERT( i < ncgroups );
    this->cgDataData()[i].ncoords = ncoords;
}

/** Close this array - call this once you have finished constructing
    the data */
void CGArrayArrayData::close()
{
    CGArrayData *cgarrays = this->cgArrayDataData();
    CGData *cgroups = this->cgDataData();
    
    //ensure that everything adds up...
    quint32 n_assigned_cgroups = 0;
    quint32 n_assigned_coords = 0;
    
    //first close out the CoordGroups
    for (quint32 i=0; i<ncgroups; ++i)
    {
        CGData &cgroup = cgroups[i];

        //tell the CoordGroup where its first point is located
        if (cgroup.nCoords() > 0)
            cgroup.coords0 = coords0 + n_assigned_coords * sizeof(Vector);
        else
            cgroup.coords0 = 0;
        
        //tell the CoordGroup where its AABox is located
        cgroup.aabox = aabox0 + i*sizeof(AABox);
        
        n_assigned_coords += cgroup.nCoords();
    }

    BOOST_ASSERT( n_assigned_coords == ncoords );

    n_assigned_coords = 0;

    //now close out the CoordGroupArrays
    for (quint32 i=0; i<ncgarrays; ++i)
    {
        CGArrayData &cgarray = cgarrays[i];
    
        if (cgarray.nCGroups() > 0)
        {
            //tell the array where its first CoordGroup is located
            cgarray.cgroup0 = cgroup0 + n_assigned_cgroups * sizeof(CoordGroup);
        
            //tell the array where its first AABox is located
            cgarray.aabox0 = aabox0 + n_assigned_cgroups * sizeof(AABox);
        
            //sum up the number of coordinates in the CoordGroups in this array
            quint32 ncgcoords = 0;
            const CoordGroup *array_groups = cgarray.constCGroupData();
        
            for (quint32 j=0; j<cgarray.nCGroups(); ++j)
            {
                ncgcoords += array_groups[j].count();
            }
        
            cgarray.ncoords = ncgcoords;
            
            if (ncgcoords > 0)
            {
                //tell the array where its first point is located
                cgarray.coords0 = coords0 + n_assigned_coords * sizeof(Vector);
                n_assigned_coords += ncgcoords;
            }
            else
                cgarray.coords0 = 0;
        }
        else
        {
            cgarray.cgroup0 = 0;
            cgarray.aabox0 = 0;
            cgarray.coords0 = 0;
        }
        
        n_assigned_cgroups += cgarray.ncgroups;
    }
    
    BOOST_ASSERT( n_assigned_coords == ncoords );
    BOOST_ASSERT( n_assigned_cgroups == ncgroups );
}

////////
//////// Implementation of CGArrayData
////////

/** Null constructor */
CGArrayData::CGArrayData()
            : this_cgarray(0), cgroup0(0), aabox0(0),
              ncgroups(0), coords0(0), ncoords(0)
{}

/** Construct, telling this array that it is at index 'this_idx' */
CGArrayData::CGArrayData(quint32 this_idx)
            : this_cgarray(this_idx), cgroup0(0), aabox0(0),
              ncgroups(0), coords0(0), ncoords(0)
{}

/** Copy constructor */
CGArrayData::CGArrayData(const CGArrayData &other)
            : this_cgarray(other.this_cgarray), 
              cgroup0(other.cgroup0), 
              aabox0(other.aabox0), ncgroups(other.ncgroups),
              coords0(other.coords0), ncoords(other.ncoords)
{}

/** Destructor */
CGArrayData::~CGArrayData()
{}

/** Return a const pointer to the start of the raw memory
    that holds all of the data of everything in the container */
const char* CGArrayData::memory() const
{
    return CGMemory::getRoot( (const char*)this, this_cgarray );
}

/** Return a pointer to the start of the raw memory
    that holds all of the data of everything in the container.
    This will detach from shared storage */
char* CGArrayData::memory()
{
    return CGMemory::getRoot( (char*)this, this_cgarray );
}

/** Increment the reference count for this object */
void CGArrayData::incref()
{
    CGMemory::incref( (char*)this, this_cgarray );
}

/** Decrement the reference count for this object - this may delete it! */
void CGArrayData::decref()
{
    CGMemory::decref( (char*)this, this_cgarray );
}

/** Detach from shared storage */
CGArrayData* CGArrayData::detach()
{
    return (CGArrayData*)( CGMemory::detach( (char*)this, this_cgarray ) );
}

/** Return a pointer to a CGArrayArrayData object that contains 
    just the data of this CoordGroupArray */
CGArrayArrayData* CGArrayData::extract() const
{
    const CGArrayArrayData *array = (const CGArrayArrayData*)(memory());
    
    if (array->nCGArrays() <= 1)
        //this is already extracted!
        return const_cast<CGArrayArrayData*>(array);
    
    char *new_storage = CGMemory::create(1, this->nCGroups(), this->nCoords());
    
    try
    {
                
    CGArrayArrayData *new_array = (CGArrayArrayData*)(new_storage);

    new_array->setNCGroupsInArray(0, this->nCGroups());
    
    for (quint32 i=0; i<this->nCGroups(); ++i)
    {
        new_array->setNPointsInCGroup(i, cGroupData()[i].count());
    }

    new_array->close();
    
    //copy the AABoxes
    AABox *output = quickCopy<AABox>( new_array->aaBoxData(),
                                      this->aaBoxData(), 
                                      this->nCGroups() );
                             
    BOOST_ASSERT( output == new_array->aaBoxData() );

    //copy the coordinates
    Vector *out2 = quickCopy<Vector>( new_array->coordsData(),
                                      this->coordsData(),
                                      this->nCoords() );
                       
    BOOST_ASSERT( out2 == new_array->coordsData() );

    return new_array;

    }
    catch(...)
    {
        delete[] new_storage;
    }
    
    return 0;
}

/** Return a pointer to the first CoordGroup in this array.
    This returns 0 if there are no CoordGroups in this array */
const CoordGroup* CGArrayData::cGroupData() const
{
    if (cgroup0 == 0)
        return 0;
    else
        return (const CoordGroup*)( memory() + cgroup0 );
}

/** Return a pointer to the AABox belonging to the first CoordGroup
    in this array. This returns 0 if there are no CoordGroups in this
    array */
const AABox* CGArrayData::aaBoxData() const
{
    if (aabox0 == 0)
        return 0;
    else
        return (const AABox*)( memory() + aabox0 );
}

/** Return a pointer to the coordinates of the first point in the 
    first CoordGroup in this array. This returns 0 if there are
    no points in this array */
const Vector* CGArrayData::coordsData() const
{
    if (coords0 == 0)
        return 0;
    else
        return (const Vector*)( memory() + coords0 );
}

/** Return a pointer to the first CoordGroup in this array.
    This returns 0 if there are no CoordGroups in this array */
const CoordGroup* CGArrayData::constCGroupData() const
{
    if (cgroup0 == 0)
        return 0;
    else
        return (const CoordGroup*)( memory() + cgroup0 );
}

/** Return a pointer to the AABox belonging to the first CoordGroup
    in this array. This returns 0 if there are no CoordGroups in this
    array */
const AABox* CGArrayData::constAABoxData() const
{
    if (aabox0 == 0)
        return 0;
    else
        return (const AABox*)( memory() + aabox0 );
}

/** Return a pointer to the coordinates of the first point in the 
    first CoordGroup in this array. This returns 0 if there are
    no points in this array */
const Vector* CGArrayData::constCoordsData() const
{
    if (coords0 == 0)
        return 0;
    else
        return (const Vector*)( memory() + coords0 );
}

/** Return a pointer to the first CoordGroup in this array.
    This returns 0 if there are no CoordGroups in this array */
CoordGroup* CGArrayData::cGroupData()
{
    if (cgroup0 == 0)
        return 0;
    else
        return (CoordGroup*)( memory() + cgroup0 );
}

/** Return a pointer to the AABox belonging to the first CoordGroup
    in this array. This returns 0 if there are no CoordGroups in this
    array */
AABox* CGArrayData::aaBoxData()
{
    if (aabox0 == 0)
        return 0;
    else
        return (AABox*)( memory() + aabox0 );
}

/** Return a pointer to the coordinates of the first point in the 
    first CoordGroup in this array. This returns 0 if there are
    no points in this array */
Vector* CGArrayData::coordsData()
{
    if (coords0 == 0)
        return 0;
    else
        return (Vector*)( memory() + coords0 );
}

/** Return the number of CoordGroups in this array */
quint32 CGArrayData::nCGroups() const
{
    return ncgroups;
}

/** Return the number of coordinates (in total) in this array */
quint32 CGArrayData::nCoords() const
{
    return ncoords;
}

////////
//////// Implementation of CGData
////////

/** Null constructor */
CGData::CGData()
       : this_cgroup(0), coords0(0), ncoords(0), aabox(0)
{}

/** Construct the data, telling it that it is at 
    index 'this_idx' in the memory storage array */
CGData::CGData(quint32 this_idx)
       : this_cgroup(this_idx), coords0(0), ncoords(0), aabox(0)
{}

/** Copy constructor */
CGData::CGData(const CGData &other)
       : this_cgroup(other.this_cgroup), coords0(other.coords0),
         ncoords(other.ncoords), aabox(other.aabox)
{}

/** Destructor */
CGData::~CGData()
{}

/** Return a pointer to a CGData that has been detached from 
    shared storage (could be a pointer to this CGData) */
CGData* CGData::detach()
{
    return (CGData*)( CGMemory::detach( (char*)this, this_cgroup ) );
}

/** Return a pointer to a CGArrayArrayData object that contains 
    just the data of this CoordGroup */
CGArrayArrayData* CGData::extract() const
{
    const CGArrayArrayData *array = (const CGArrayArrayData*)(memory());
    
    if (array->nCGroups() <= 1)
        //this is already extracted!
        return const_cast<CGArrayArrayData*>(array);
    
    char *new_storage = CGMemory::create(1, 1, this->nCoords());
    
    try
    {
                
    //ok, we need to extract!
    CGArrayArrayData *new_array = (CGArrayArrayData*)(new_storage);

    new_array->setNCGroupsInArray(0, 1);
    new_array->setNPointsInCGroup(0, this->nCoords());
    new_array->close();
    
    //copy the AABox
    new_array->aaBoxData()[0] = *(this->aaBox());

    //copy the coordinates
    Vector *output = quickCopy<Vector>( new_array->coordsData(),
                                        this->coordsData(),
                                        this->nCoords() );
                       
    BOOST_ASSERT( output == new_array->coordsData() );

    return new_array;

    }
    catch(...)
    {
        delete[] new_storage;
    }
    
    return 0;
}

/** Return a const pointer to the start of the storage
    array that contains this CGData */
const char* CGData::memory() const
{
    return CGMemory::getRoot( (const char*)this, this_cgroup );
}

/** Return a pointer to the start of the storage array that
    contains this CGData - this will detach from shared storage */
char* CGData::memory()
{
    return CGMemory::getRoot( (char*)this, this_cgroup );
}

/** Increment the reference count for this object */
void CGData::incref()
{
    CGMemory::incref( (char*)this, this_cgroup );
}

/** Decrement the reference count for this object - this may delete it! */
void CGData::decref()
{
    CGMemory::decref( (char*)this, this_cgroup );
}

/** Return a pointer to the coordinates of the first point in this
    CoordGroup. This returns 0 if there are no points in this CoordGroup */
const Vector* CGData::coordsData() const
{
    if (coords0 == 0)
        return 0;
    else
        return (const Vector*)( memory() + coords0 );
}

/** Return a pointer to the AABox that surrounds all of the points
    in this CoordGroup */
const AABox* CGData::aaBox() const
{
    return (const AABox*)( memory() + aabox );
}

/** Return a pointer to the coordinates of the first point in this
    CoordGroup. This returns 0 if there are no points in this CoordGroup */
Vector* CGData::coordsData()
{
    if (coords0 == 0)
        return 0;
    else
        return (Vector*)( memory() + coords0 );
}

/** Return a pointer to the AABox that surrounds all of the points
    in this CoordGroup */
AABox* CGData::aaBox()
{
    return (AABox*)( memory() + aabox );
}

/** Return the number of points in this CoordGroup */
quint32 CGData::nCoords() const
{
    return ncoords;
}

////////
//////// Implementation of CoordGroupBase
////////

static CGSharedPtr<CGArrayArrayData> shared_null;

static const CGSharedPtr<CGArrayArrayData>& getSharedNull()
{
    if (shared_null.constData() == 0)
    {
        CGSharedPtr<CGArrayArrayData> my_null = (CGArrayArrayData*)( CGMemory::create(0,0,0) );

        auto mutex = SireBase::detail::get_shared_null_mutex();
        tbb::spin_mutex::scoped_lock lock(*mutex);
    
        if (shared_null.constData() == 0)
        {
            shared_null = my_null;
        }
        else
            lock.release();
    }
    
    return shared_null;
}

static CGSharedPtr<CGData> getSharedNullCoordGroup()
{
    const CGSharedPtr<CGArrayArrayData> &array = ::getSharedNull();
    
    return CGSharedPtr<CGData>( array->nullCGroup() );
}

/** Null constructor */
CoordGroupBase::CoordGroupBase()
               : d( ::getSharedNullCoordGroup() )
{}

/** Private constructor used only by CGData */
CoordGroupBase::CoordGroupBase(CGData *data)
{
    d.weakAssign(data);
}

static CGSharedPtr<CGData> createCoordGroup(quint32 size)
{
    if (size == 0)
        return ::getSharedNullCoordGroup();

    //construct space for 1 array of 1 CoordGroup of size atoms
    char *storage = CGMemory::create(1, 1, size);
        
    CGArrayArrayData *array = (CGArrayArrayData*)storage;
        
    array->setNCGroupsInArray(0, 1);
    array->setNPointsInCGroup(0, size);
        
    //close the array - this validates that it is complete
    array->close();
    
    return CGSharedPtr<CGData>( &(array->cgDataData()[0]) );
}

/** Construct a CoordGroup that can hold up to 'size'
    coordinates, each of value 'value' */
CoordGroupBase::CoordGroupBase(quint32 size, const Vector &value)
                : d( createCoordGroup(size) )
{
    if (size > 0 and not value.isZero())
    {
        Vector *coords = d->coordsData();
        
        //initialise all of the coordinates to the specified value
        for (quint32 i=0; i<size; ++i)
        {
            coords[i] = value;
        }
        
        //update the AABox
        *(d->aaBox()) = AABox(value);
    }
}

/** Construct a CoordGroup that holds the 'size' points
    copied from the array 'values' */
CoordGroupBase::CoordGroupBase(quint32 size, const Vector *values)
                : d( createCoordGroup(size) )
{
    if (size > 0)
    {
        Vector *coords = d->coordsData();
        Vector *output = quickCopy<Vector>(coords, values, size);

        BOOST_ASSERT( output == coords );
        
        *(d->aaBox()) = AABox(*this);
    }
}

/** Construct a CoordGroup from all of the coordinates in the 
    passed CoordGroupArray */
CoordGroupBase::CoordGroupBase(const CoordGroupArray &cgarray)
               : d( getSharedNullCoordGroup() )
{
    if (cgarray.nCoordGroups() == 1)
    {
        this->operator=(cgarray.constData()[0]);
    }
    else if (cgarray.nCoords() > 0)
    {
        d = createCoordGroup(cgarray.nCoords());
    
        Vector *coords = d->coordsData();
        Vector *output = quickCopy<Vector>(coords, cgarray.constCoordsData(), 
                                           cgarray.nCoords());
        
        BOOST_ASSERT( output == coords );
        
        *(d->aaBox()) = AABox(cgarray);
    }
}

/** Construct a CoordGroup from all of the coordinates in all
    of the arrays in the passed set of CoordGroupArrays */
CoordGroupBase::CoordGroupBase(const CoordGroupArrayArray &cgarrays)
               : d( getSharedNullCoordGroup() )
{
    if (cgarrays.nCoordGroups() == 1)
    {
        this->operator=(cgarrays.constCoordGroupData()[0]);
    }
    else if (cgarrays.nCoords() > 0)
    {
        d = createCoordGroup(cgarrays.nCoords());
    
        Vector *coords = d->coordsData();
        Vector *output = quickCopy<Vector>(coords, cgarrays.constCoordsData(),
                                           cgarrays.nCoords());

        BOOST_ASSERT( output == coords );
                 
        *(d->aaBox()) = AABox(cgarrays);
    }
}

/** Construct a CoordGroup from the specified coordinates */
CoordGroupBase::CoordGroupBase(const QVector<Vector> &coordinates)
                : d( createCoordGroup(coordinates.count()) )
{
    //how many coordinates are there?
    quint32 ncoords = coordinates.count();

    if (ncoords > 0)
    {
        Vector *coords = d->coordsData();
        
        //block copy the coordinates
        Vector *output = quickCopy<Vector>( coords, coordinates.constData(),
                                            ncoords );

        BOOST_ASSERT( output == coords );
        
        //now update the AABox
        *(d->aaBox()) = AABox(coordinates);
    }
}

/** Copy constructor - this is fast as CoordGroupBase is implicitly
    shared */
CoordGroupBase::CoordGroupBase(const CoordGroupBase &other)
               : d(other.d)
{}

/** Destructor - this will only delete the data if this
    is the only reference to it */
CoordGroupBase::~CoordGroupBase()
{}

/** Copy assignment operator - CoordGroupBase is implicitly
    shared so this will be fast */
CoordGroupBase& CoordGroupBase::operator=(const CoordGroupBase &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool CoordGroupBase::operator==(const CoordGroupBase &other) const
{
    if ( d.constData() == other.d.constData() )
    {
        return true;
    }
    else if (this->count() != other.count())
    {
        return false;
    }
    else
    {
        //compare the coordinates
        const Vector *my_coords = constData();
        const Vector *other_coords = other.constData();

        int sz = size();
        for (int i=0; i<sz; ++i)
        {
            if (my_coords[i] != other_coords[i])
                return false;
        }

        return true;
    }
}

/** Comparison operator */
bool CoordGroupBase::operator!=(const CoordGroupBase &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this CoordGroup */
QString CoordGroupBase::toString() const
{
    if (this->isEmpty())
        return QObject::tr("CoordGroup::null");
        
    QStringList lines;
    
    for (int i=0; i<this->count(); ++i)
    {
        const Vector &c = this->constData()[i];
        lines.append( QString("[%1, %2, %3]").arg(c.x()).arg(c.y()).arg(c.z()) );
    }
    
    return QObject::tr("CoordGroup( aaBox() == %1,\n    %2\n         )")
                .arg(this->aaBox().toString(), lines.join("\n    "));
}

/** Return whether 'other' may be different to this group - this uses
    a simple comparison of the memory addresses of the storage of
    these two groups to see if they are definitely the same, or maybe
    different */
bool CoordGroupBase::maybeDifferent(const CoordGroupBase &other) const
{
    return d.constData() != other.d.constData();
}

/** Return whether this group is empty (has no coordinates) */
bool CoordGroupBase::isEmpty() const
{
    return d->nCoords() == 0;
}

/** Return the enclosing AABox */
const AABox& CoordGroupBase::aaBox() const
{
    return *(d->aaBox());
}

/** Return a raw pointer to the array containing all
    of the coordinates */
const Vector* CoordGroupBase::constData() const
{
    return d->coordsData();
}

/** Return a raw pointer to the array containing all
    of the coordinates */
const Vector* CoordGroupBase::data() const
{
    return d->coordsData();
}

/** Return the number of coordinates in this group */
int CoordGroupBase::count() const
{
    return d->nCoords();
}

/** Return the number of coordinates in this group */
int CoordGroupBase::size() const
{
    return d->nCoords();
}

/** Return an array containing just the coordinates in this group */
QVector<Vector> CoordGroupBase::toVector() const
{
    const int ncoords = d->nCoords();
    
    if (ncoords == 0)
        return QVector<Vector>();
        
    else
    {
        QVector<Vector> coords(ncoords);
        
        quickCopy<Vector>(coords.data(), this->constData(), ncoords);
        
        return coords;
    }
}

static SIRE_ALWAYS_INLINE void assertSame(quint32 sz, quint32 other)
{
    if (sz != other)
        throw SireError::incompatible_error( QObject::tr(
            "The number of coordinates in this group is %1, while "
            "the number of coordinates in the other is %2!")
                .arg(sz).arg(other), CODELOC );
}

/** Assert that the array of coordinates contains the same number of
    coordinates as are in this group

    \throw SireError::incompatible_error
*/
void CoordGroupBase::assertSameSize(const QVector<Vector> &coordinates) const
{
    assertSame(this->count(), coordinates.count());
}

/** Assert that the CoordGroup 'other' contains the same number of
    coordinates as are in this group

    \throw SireError::incompatible_error
*/
void CoordGroupBase::assertSameSize(const CoordGroupBase &other) const
{
    assertSame(this->count(), other.count());
}

/** Assert that 'i' is a valid index in this group

    \throw SireError::invalid_index
*/
void CoordGroupBase::assertValidIndex(quint32 i) const
{
    if (i >= quint32(this->size()))
        throw SireError::invalid_index( QObject::tr(
            "Invalid index: %1 - maximum index == %2.")
                .arg(i).arg(this->size()), CODELOC );
}

/** Return the 'ith' coordinate of this group

    \throw SireError::invalid_index
*/
const Vector& CoordGroupBase::at(quint32 i) const
{
    assertValidIndex(i);
    return d->coordsData()[i];
}

/** Return the 'ith' coordinate of this group

    \throw SireError::invalid_index
*/
const Vector& CoordGroupBase::operator[](quint32 i) const
{
    assertValidIndex(i);
    return d->coordsData()[i];
}

////////
//////// Implementation of CoordGroup
////////

static const RegisterMetaType<CoordGroup> r_cgroup(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const CoordGroup &cgroup)
{
    writeHeader(ds, r_cgroup, 1);

    SharedDataStream sds(ds);

    sds << CoordGroupArrayArray(cgroup);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       CoordGroup &cgroup)
{
    VersionID v = readHeader(ds, r_cgroup);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
    
        CoordGroupArrayArray array;
        sds >> array;

        BOOST_ASSERT( array.coordGroupData() != 0 );
    
        cgroup = array.coordGroupData()[0];
    }
    else
        throw version_error(v, "1", r_cgroup, CODELOC);
    
    return ds;
}

/** Null constructor */
CoordGroup::CoordGroup() : CoordGroupBase()
{}

/** Private constructor used only by CGData */
CoordGroup::CoordGroup(CGData *data) : CoordGroupBase(data)
{}

/** Construct a CoordGroup that holds 'size' coordinates */
CoordGroup::CoordGroup(quint32 size)
           : CoordGroupBase(size)
{}

/** Construct a CoordGroup that holds 'size' coordinates,
    all of which have the value 'value' */
CoordGroup::CoordGroup(quint32 size, const Vector &value)
           : CoordGroupBase(size,value)
{}

/** Construct a CoordGroup that holds the 'size' points
    copied from the array 'values' */
CoordGroup::CoordGroup(quint32 size, const Vector *values)
           : CoordGroupBase(size, values)
{}

/** Construct a CoordGroup from all of the coordinates contained
    in the passed CoordGroup array */
CoordGroup::CoordGroup(const CoordGroupArray &cgarray)
           : CoordGroupBase(cgarray)
{}

/** Construct a CoordGroup from all of the coordinates in the
    passed collection of CoordGroup arrays */
CoordGroup::CoordGroup(const CoordGroupArrayArray &cgarrays)
           : CoordGroupBase(cgarrays)
{}

/** Construct a CoordGroup from the array of passed coordinates */
CoordGroup::CoordGroup(const QVector<Vector> &points)
           : CoordGroupBase(points)
{}

/** Construct from a CoordGroupEditor (only CoordGroupEditor::commit()
    calls this function, and in this case it has already made sure that
    the AABox is up-to-date) */
CoordGroup::CoordGroup(const CoordGroupEditor &editor)
            : CoordGroupBase(editor)
{}

/** Copy constructor */
CoordGroup::CoordGroup(const CoordGroup &other)
           : CoordGroupBase(other)
{}

/** Destructor */
CoordGroup::~CoordGroup()
{}

void CoordGroup::throwInvalidCountError(uint nats0, uint nats1)
{
    throw SireError::invalid_state( QObject::tr(
        "You can only split a CoordGroup into groups if the total number "
        "of atoms in the single CoordGroup and the group of CoordGroups "
        "is the same. It isn't in this case! %1 vs. %2.")
            .arg(nats0).arg(nats1), CODELOC );
}

/** Copy assignment operator */
CoordGroup& CoordGroup::operator=(const CoordGroup &other)
{
    CoordGroupBase::operator=(other);
    return *this;
}

/** Copy assignment from a CoordGroupEditor */
CoordGroup& CoordGroup::operator=(CoordGroupEditor &other)
{
    return this->operator=( other.commit() );
}

/** Return an editor that can be used to edit the
    coordinates in this group */
CoordGroupEditor CoordGroup::edit() const
{
    return CoordGroupEditor(*this);
}

const char* CoordGroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CoordGroup>() );
}

////////
//////// Implementation of CoordGroupEditor
////////

/** Null constructor */
CoordGroupEditor::CoordGroupEditor() 
                  : CoordGroupBase(), needsupdate(false)
{}

/** Construct from a CoordGroup */
CoordGroupEditor::CoordGroupEditor(const CoordGroup &other)
                 : CoordGroupBase(other), needsupdate(false)
{
    if (other.count() != 0)
        d = other.d->extract()->cGroupData()[0].d;
}

/** Copy constructor */
CoordGroupEditor::CoordGroupEditor(const CoordGroupEditor &other)
                  : CoordGroupBase(other), needsupdate(other.needsupdate)
{}

/** Destructor */
CoordGroupEditor::~CoordGroupEditor()
{}

/** Copy assignment operator */
CoordGroupEditor& CoordGroupEditor::operator=(const CoordGroupEditor &other)
{
    CoordGroupBase::operator=(other);
    needsupdate = other.needsupdate;
    
    return *this;
}

/** Assignment operator */
CoordGroupEditor& CoordGroupEditor::operator=(const CoordGroup &other)
{
    if (other.count() > 0)
        d = other.d->extract()->cGroupData()[0].d;
    else
        d = 0;
    
    needsupdate = false;
    
    return *this;
}

/** Return a modifiable reference to the 'ith' coordinate in the group
    - this will throw an exception if 'i' refers to an invalid index

    \throw SireError::invalid_index
*/
Vector& CoordGroupEditor::operator[](quint32 i)
{
    CoordGroupBase::assertValidIndex(i);
    needsupdate = true;

    return d->coordsData()[i];
}

/** Return a modifiable pointer to the array of coordinates -
    *do not delete this pointer* */
Vector* CoordGroupEditor::data()
{
    needsupdate = true;
    return d->coordsData();
}

/** Translate this CoordGroup by 'delta' */
CoordGroupEditor& CoordGroupEditor::translate(const Vector &delta)
{
    if (not delta.isZero())
    {
        //translate all of the coordinates
        quint32 ncoords = this->count();
        
        Vector *coords = d->coordsData();

        for (quint32 i=0; i<ncoords; ++i)
        {
            coords[i] += delta;
        }

        //now translate the AABox
        if (not needsupdate)
        {
            AABox new_box = this->aaBox();
            new_box.translate(delta);
            *(d->aaBox()) = new_box;
        }
    }

    return *this;
}

/** Translate the 'ith' point in the group by 'delta'

    \throw SireError::invalid_index
*/
CoordGroupEditor& CoordGroupEditor::translate(quint32 i, const Vector &delta)
{
    if (not delta.isZero())
        this->operator[](i) += delta;

    needsupdate = true;

    return *this;
}

/** Rotate (and perhaps scale and shear) this group by the matrix 'rotmat'
    about the point 'point' */
CoordGroupEditor& CoordGroupEditor::rotate(const Matrix &rotmat, const Vector &point)
{
    quint32 sz = this->count();
    Vector *coords = d->coordsData();

    for (quint32 i=0; i<sz; ++i)
    {
        coords[i] = SireMaths::rotate( coords[i], rotmat, point );
    }

    needsupdate = true;

    return *this;
}

/** Transform this group by the transformation 't' */
CoordGroupEditor& CoordGroupEditor::transform(const Transform &t)
{
    quint32 sz = this->count();
    Vector *coords = d->coordsData();

    coords = t.apply(coords, sz);

    needsupdate = true;

    return *this;
}

/** Rotate the 'ith' point in the group using the matrix 'rotmat' about the
    point 'point'

    \throw SireError::index
*/
CoordGroupEditor& CoordGroupEditor::rotate(quint32 i, const Matrix &rotmat, 
                                           const Vector &point)
{
    Vector &coord = this->operator[](i);
    coord = SireMaths::rotate(coord, rotmat, point);
    needsupdate = true;

    return *this;
}

/** Rotate the 'ith' point in the group using the matrix 'rotmat' about the
    point 'point'

    \throw SireError::index
*/
CoordGroupEditor& CoordGroupEditor::transform(quint32 i, const Transform &t)
{
    Vector &coord = this->operator[](i);
    coord = t.apply(coord);
    needsupdate = true;

    return *this;
}

/** Rotate this group by the Quaternion 'quat' about the point 'point' */
CoordGroupEditor& CoordGroupEditor::rotate(const Quaternion &quat, 
                                             const Vector &point)
{
    return this->rotate(quat.toMatrix(), point);
}

/** Rotate the 'ith' point in the group using the quaternion 'quat' about the
    point 'point'

    \throw SireError::index
*/
CoordGroupEditor& CoordGroupEditor::rotate(quint32 i, const Quaternion &quat, 
                                           const Vector &point)
{
    return this->rotate(i, quat.toMatrix(), point);
}

/** Set the coordinates of the CoordGroup to 'newcoords' - this must
    have the same number of points as this CoordGroup or an exception
    will be thrown

    \throw SireError::incompatible_error
*/
CoordGroupEditor& CoordGroupEditor::setCoordinates(const QVector<Vector> &newcoords)
{
    CoordGroupBase::assertSameSize(newcoords);

    int sz = this->count();
    Vector *coords = d->coordsData();
    const Vector *newcoords_array  = newcoords.constData();

    Vector *output = quickCopy<Vector>(coords, newcoords_array, sz);
    needsupdate = true;

    BOOST_ASSERT(output == coords);

    return *this;
}

/** Set the coordinates of the ith atom to 'newcoords'

    \throw SireError::invalid_index
*/
CoordGroupEditor& CoordGroupEditor::setCoordinates(quint32 i, const Vector &newcoords)
{
    if (i >= quint32(this->count()))
        throw SireError::invalid_index( QObject::tr(
            "Cannot update the coordinates of point %1 as the number "
            "of points available is just %2.")
                .arg(i).arg(this->count()), CODELOC );
    
    Vector *coords = d->coordsData();
    coords[i] = newcoords;
    
    needsupdate = true;
    
    return *this;
}

/** Set the coordinates of the CoordGroup to 'newcoords' - this
    must have the same number of points as this CoordGroup or an
    exception will be thrown

    \throw SireError::incompatible_error
*/
CoordGroupEditor& CoordGroupEditor::setCoordinates(const CoordGroupBase &newcoords)
{
    if ( this->constData() != newcoords.constData() )
    {
        CoordGroupBase::assertSameSize(newcoords);

        int sz = this->count();
        Vector *coords = d->coordsData();
        const Vector *newcoords_array  = newcoords.constData();

        Vector *output = quickCopy<Vector>(coords, newcoords_array, sz);
        needsupdate = true;

        BOOST_ASSERT(output == coords);
    }

    return *this;
}

/** Map all of these coordinates into the axis set 'axes' */
CoordGroupEditor& CoordGroupEditor::mapInto(const AxisSet &axes)
{
    int sz = this->count();
    
    Vector *coords = d->coordsData();
    
    for (int i=0; i<sz; ++i)
    {
        coords[i] = axes.fromIdentity(coords[i]);
    }
    
    needsupdate = true;
    
    return *this;
}

/** Map the coordinates of the ith point into the axis set 'axes'

    \throw SireError::invalid_index
*/
CoordGroupEditor& CoordGroupEditor::mapInto(quint32 i, const AxisSet &axes)
{
    if (i >= quint32(this->count()))
        throw SireError::invalid_index( QObject::tr(
            "Cannot update the coordinates of point %1 as the number "
            "of points available is just %2.")
                .arg(i).arg(this->count()), CODELOC );
    
    Vector *coords = d->coordsData();
    coords[i] = axes.fromIdentity(coords[i]);
    
    needsupdate = true;
    
    return *this;
}

/** Change the coordinate frame of all of the coordinates from
    'from_frame' to 'to_frame'. This has the same effect
    as 'mapInto' if 'from_frame' is the unit cartesian set */
CoordGroupEditor& CoordGroupEditor::changeFrame(const AxisSet &from_frame,
                                                const AxisSet &to_frame)
{
    int sz = this->count();
    
    Vector *coords = d->coordsData();

    //precalculate what we can...
    Matrix mat = from_frame.invMatrix() * to_frame.matrix();
    
    for (int i=0; i<sz; ++i)
    {
        coords[i] = (mat * (coords[i] - from_frame.origin())) + to_frame.origin();
    }
    
    needsupdate = true;
    
    return *this;
}

/** Map the coordinates of the ith point from the frame
    'from_frame' to the frame 'to_frame'

    \throw SireError::invalid_index
*/
CoordGroupEditor& CoordGroupEditor::changeFrame(quint32 i,
                                                const AxisSet &from_frame,
                                                const AxisSet &to_frame)
{
    if (i >= quint32(this->count()))
        throw SireError::invalid_index( QObject::tr(
            "Cannot update the coordinates of point %1 as the number "
            "of points available is just %2.")
                .arg(i).arg(this->count()), CODELOC );
    
    Vector *coords = d->coordsData();
    coords[i] = from_frame.toFrame(to_frame,coords[i]);
    
    needsupdate = true;
    
    return *this;
}

/** Return a CoordGroup which is a copy of this group. */
CoordGroup CoordGroupEditor::commit()
{
    //return a copy of this CoordGroup
    if (needsupdate)
    {
        *(d->aaBox()) = AABox(*this);
        needsupdate = false;
    }
    
    return CoordGroup(*this);
}

const char* CoordGroupEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CoordGroupEditor>() );
}

////////
//////// Implementation of CoordGroupArray
////////

static const RegisterMetaType<CoordGroupArray> r_cgarray(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const CoordGroupArray &cgarray)
{
    writeHeader(ds, r_cgarray, 1);
    
    SharedDataStream sds(ds);
    
    sds << CoordGroupArrayArray(cgarray);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       CoordGroupArray &cgarray)
{
    VersionID v = readHeader(ds, r_cgarray);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        CoordGroupArrayArray array;
        sds >> array;

        BOOST_ASSERT( array.data() != 0 );
        
        cgarray = array.data()[0];
    }
    else
        throw version_error(v, "1", r_cgarray, CODELOC);
        
    return ds;
}

static CGSharedPtr<CGArrayData> getSharedNullCGArray()
{
    const CGSharedPtr<CGArrayArrayData> &array = ::getSharedNull();
    
    return CGSharedPtr<CGArrayData>( array->nullArray() );
}

/** Null constructor */
CoordGroupArray::CoordGroupArray()
                : d( ::getSharedNullCGArray() )
{}

/** Private constructor used only by CGArrayData */
CoordGroupArray::CoordGroupArray(CGArrayData *data)
{
    d.weakAssign(data);
}

/** Construct an array that holds just the passed CoordGroup */
CoordGroupArray::CoordGroupArray(const CoordGroup &cgroup)
{
    d = cgroup.d->extract()->cgArrayData()[0].d;
}

static CGSharedPtr<CGArrayData> createCGArray(quint32 ngroups, quint32 ncoords)
{
    if (ngroups == 0)
        return ::getSharedNullCGArray();

    //construct space for 1 array of ngroups CoordGroups of ncoords coordinates
    char *storage = CGMemory::create(1, ngroups, ncoords);
        
    CGArrayArrayData *array = (CGArrayArrayData*)storage;
    
    return CGSharedPtr<CGArrayData>( &(array->cgArrayDataData()[0]) );
}

/** Construct from an array of CoordGroups */
CoordGroupArray::CoordGroupArray(const QVector<CoordGroup> &cgroups)
{
    //count the number of groups and coordinates
    quint32 ncoords = 0;
    
    const CoordGroup *cgroups_array = cgroups.constData();
    
    for (int i=0; i<cgroups.count(); ++i)
    {
        ncoords += cgroups_array[i].count();
    }
    
    //now create space for the CoordGroups
    d = ::createCGArray(cgroups.count(), ncoords);
    
    if (cgroups.isEmpty())
        return;
    
    CGArrayArrayData *cgarrayarray = (CGArrayArrayData*)(d->memory());
    
    cgarrayarray->setNCGroupsInArray(0, cgroups.count());
    
    //tell the array about the sizes of the groups
    for (int i=0; i<cgroups.count(); ++i)
    {
        const CoordGroup &cgroup = cgroups_array[i];
        cgarrayarray->setNPointsInCGroup(i, cgroup.count());
    }
    
    cgarrayarray->close();
    
    Vector *coords = cgarrayarray->coordsData();
    AABox *aabox = cgarrayarray->aaBoxData();
    
    //finally, copy the data into the arrays
    for (int i=0; i<cgroups.count(); ++i)
    {
        const CoordGroup &cgroup = cgroups_array[i];
        
        //copy the coordinates
        Vector *output = quickCopy<Vector>(coords, cgroup.constData(),
                                           cgroup.count());
                                
        BOOST_ASSERT(output == coords);
        
        coords += cgroup.count();
        
        //now copy the AABox
        aabox[i] = cgroup.aaBox();
    }
}

/** Construct from a double-vector */
CoordGroupArray::CoordGroupArray(const QVector< QVector<Vector> > &points)
{
    if (points.isEmpty())
        return;
        
    QVector<CoordGroup> cgroups(points.count());
    
    for (int i=0; i<points.count(); ++i)
    {
        cgroups[i] = CoordGroup(points[i]);
    }
    
    this->operator=( CoordGroupArray(cgroups) );
}

/** Construct from a pair of CoordGroupArrays */
CoordGroupArray::CoordGroupArray(const CoordGroupArray &array0,
                                 const CoordGroupArray &array1)
{
    if (array0.isEmpty())
        d = array1.d;
    else if (array1.isEmpty())
        d = array0.d;
    else
    {
        //now create space for the CoordGroups
        d = ::createCGArray(array0.nCoordGroups() + array1.nCoordGroups(),
                            array0.nCoords() + array1.nCoords());
    
        CGArrayArrayData *cgarrayarray = (CGArrayArrayData*)(d->memory());
    
        cgarrayarray->setNCGroupsInArray(0, array0.nCoordGroups() +
                                            array1.nCoordGroups());
    
        for (int i=0; i<array0.nCoordGroups(); ++i)
        {
            cgarrayarray->setNPointsInCGroup(i, array0.constData()[i].count());
        }
        
        for (int i=0; i<array1.nCoordGroups(); ++i)
        {
            cgarrayarray->setNPointsInCGroup(array0.nCoordGroups() + i,
                                             array1.constData()[i].count());
        }
    
        cgarrayarray->close();
    
        Vector *coords = cgarrayarray->coordsData();
        AABox *aabox = cgarrayarray->aaBoxData();
        
        //copy the coordinates
        {
            Vector *output = quickCopy<Vector>(coords, array0.constCoordsData(),
                                               array0.nCoords());
                                
            BOOST_ASSERT(output == coords);
        
            coords += array0.nCoords();
        
            output = quickCopy<Vector>(coords, array1.constCoordsData(),
                                       array1.nCoords());
        
            BOOST_ASSERT(output == coords);
        }
        
        //copy the AABoxes
        {
            AABox *output = quickCopy<AABox>(aabox, array0.constAABoxData(),
                                             array0.nCoordGroups());
                                             
            BOOST_ASSERT(output == aabox);
        
            aabox += array0.nCoordGroups();
            
            output = quickCopy<AABox>(aabox, array1.constAABoxData(),
                                      array1.nCoordGroups());
                                      
            BOOST_ASSERT(output == aabox);
        }
    }
}

/** Copy constructor */
CoordGroupArray::CoordGroupArray(const CoordGroupArray &other)
                : d(other.d)
{}

/** Destructor */
CoordGroupArray::~CoordGroupArray()
{}

/** Copy assignment operator */
CoordGroupArray& CoordGroupArray::operator=(const CoordGroupArray &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool CoordGroupArray::operator==(const CoordGroupArray &other) const
{
    if (d.constData() == other.d.constData())
        return true;
        
    for (quint32 i=0; i<d->nCGroups(); ++i)
    {
        if (this->data()[i] != other.data()[i])
            return false;
    }
    
    return true;
}

/** Comparison operator */
bool CoordGroupArray::operator!=(const CoordGroupArray &other) const
{
    return not this->operator==(other);
}

/** Merge this array of CoordGroups back into a single CoordGroup */
CoordGroup CoordGroupArray::merge() const
{
    return CoordGroup(*this);
}

/** Assert that the index i points to a valid CoordGroup 

    \throw SireError::invalid_index
*/
void CoordGroupArray::assertValidCoordGroup(quint32 i) const
{
    if (i >= d->nCGroups())
        throw SireError::invalid_index( QObject::tr(
            "Cannot access the CoordGroup at index %1. The number "
            "of CoordGroups in this array is just %2.")
                .arg(i).arg(d->nCGroups()), CODELOC );
}

/** Assert that the index i points to a valid CoordGroup 

    \throw SireError::invalid_index
*/
void CoordGroupArray::assertValidIndex(quint32 i) const
{
    this->assertValidCoordGroup(i);
}

/** Assert that the index i points a valid coordinate

    \throw SireError::invalid_index
*/
void CoordGroupArray::assertValidCoordinate(quint32 i) const
{
    if (i >= d->nCoords())
        throw SireError::invalid_index( QObject::tr(
            "Cannot access the coordinate at index %1. The number "
            "of coordinates in this array is just %2.")
                .arg(i).arg(d->nCoords()), CODELOC );
}

/** Return a reference to the ith CoordGroup

    \throw SireError::invalid_index
*/
const CoordGroup& CoordGroupArray::operator[](quint32 i) const
{
    assertValidCoordGroup(i);
    
    return d->cGroupData()[i];
}

/** Return a reference to the ith CoordGroup

    \throw SireError::invalid_index
*/
const CoordGroup& CoordGroupArray::at(quint32 i) const
{
    return this->operator[](i);
}

/** Return whether or not this array is empty */
bool CoordGroupArray::isEmpty() const
{
    return d->nCGroups() == 0;
}

/** Return a string representation of this array */
QString CoordGroupArray::toString() const
{
    if (this->isEmpty())
        return QObject::tr("CoordGroupArray::null");
        
    QStringList lines;
    
    for (int i=0; i<this->count(); ++i)
    {
        lines.append( this->constData()[i].toString() );
    }
    
    return QObject::tr("CoordGroupArray( aaBox() == %1,\n%2\n)")
                .arg(this->aaBox().toString(), lines.join("\n"));
}

/** Return the number of CoordGroups in this array */
int CoordGroupArray::count() const
{
    return d->nCGroups();
}

/** Return the number of CoordGroups in this array */
int CoordGroupArray::size() const
{
    return d->nCGroups();
}

/** Return the number of CoordGroups in this array */
int CoordGroupArray::nCoordGroups() const
{
    return d->nCGroups();
}

/** Return the number of coordinates in this array */
int CoordGroupArray::nCoords() const
{
    return d->nCoords();
}

/** Return an AABox that complete encompasses all of the CoordGroups 
    in this array */
AABox CoordGroupArray::aaBox() const
{
    int ncgroups = this->count();

    if (ncgroups == 1)
        return this->constData()[0].aaBox();
    else if (ncgroups == 0)
        return AABox();
    else
    {
        const CoordGroup *coords_array = this->constData();
        
        Vector mincoords = coords_array[0].aaBox().minCoords();
        Vector maxcoords = coords_array[0].aaBox().maxCoords();
        
        for (int i=1; i<ncgroups; ++i)
        {
            mincoords.setMin(coords_array[i].aaBox().minCoords());
            maxcoords.setMax(coords_array[i].aaBox().maxCoords());
        }
        
        return AABox::from(mincoords, maxcoords);
    }
}

/** Return a raw pointer to the array of CoordGroups */
const CoordGroup* CoordGroupArray::data() const
{
    return d->cGroupData();
}

/** Return a raw pointer to the array of CoordGroups */
const CoordGroup* CoordGroupArray::constData() const
{
    return d->cGroupData();
}

/** Return a raw pointer to the array of all coordinates of all
    of the CoordGroups in this array */
const Vector* CoordGroupArray::coordsData() const
{
    return d->coordsData();
}

/** Return a raw pointer to the array of all coordinates of all
    of the CoordGroups in this array */
const Vector* CoordGroupArray::constCoordsData() const
{
    return d->coordsData();
}

/** Return a raw pointer to the array of all AABoxes for
    all of the CoordGroups in this array */
const AABox* CoordGroupArray::aaBoxData() const
{
    return d->aaBoxData();
}

/** Return a raw pointer to the array of all AABoxes for
    all of the CoordGroups in this array */
const AABox* CoordGroupArray::constAABoxData() const
{
    return d->aaBoxData();
}

/** Append the passed CoordGroup onto the end of this array */
void CoordGroupArray::append(const CoordGroup &group)
{
    this->operator=( CoordGroupArray(*this, group) );
}

/** Append the passed CoordGroups onto the end of this array */
void CoordGroupArray::append(const CoordGroupArray &groups)
{
    this->operator=( CoordGroupArray(*this, groups) );
}

/** Internal function used to remove the CoordGroups that are 'true' in
    the passed list */
void CoordGroupArray::pvt_remove(const QVarLengthArray<bool> &to_remove)
{
    BOOST_ASSERT( this->nCoordGroups() == to_remove.count() );
    
    const CoordGroup *cgroups_array = this->constData();
    
    //count the number of groups and coordinates
    quint32 ncoords = 0;
    quint32 ngroups = 0;
    
    for (int i=0; i<to_remove.count(); ++i)
    {
        if (not to_remove[i])
        {
            ncoords += cgroups_array[i].count();
            ++ngroups;
        }
    }
    
    if (ngroups == 0)
    {
        this->operator=( CoordGroupArray() );
        return;
    }
    
    //now create space for the CoordGroups
    CGSharedPtr<detail::CGArrayData> d2 = ::createCGArray(ngroups, ncoords);
    
    CGArrayArrayData *cgarrayarray = (CGArrayArrayData*)(d2->memory());
    
    cgarrayarray->setNCGroupsInArray(0, ngroups);
    
    int igroup = 0;
    
    for (int i=0; i<to_remove.count(); ++i)
    {
        if (not to_remove[i])
        {
            cgarrayarray->setNPointsInCGroup(igroup, cgroups_array[i].count());
            ++igroup;
        }
    }
    
    cgarrayarray->close();
    
    Vector *coords = cgarrayarray->coordsData();
    AABox *aabox = cgarrayarray->aaBoxData();
    
    //finally, copy the data into the arrays
    for (int i=0; i<to_remove.count(); ++i)
    {
        if (not to_remove[i])
        {
            const CoordGroup &cgroup = cgroups_array[i];
        
            //copy the coordinates
            Vector *output = quickCopy<Vector>(coords, cgroup.constData(),
                                               cgroup.count());
                                
            BOOST_ASSERT(output == coords);
        
            coords += cgroup.count();
        
            //now copy the AABox
            *aabox = cgroup.aaBox();
            ++aabox;
        }
    }
    
    d = d2;
}

/** Remove 'count' CoordGroups from the array, starting with the ith CoordGroup */
void CoordGroupArray::remove(quint32 start, int count)
{
    if (start == 0 and count == this->nCoordGroups())
    {
        this->operator=( CoordGroupArray() );
        return;
    }

    if (count < 0 or int(start)+count > this->nCoordGroups())
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot remove elements %1 to %2 as this is an invalid index "
                "range for a CoordGroupArray of size %3.")
                    .arg(start).arg(start+count)
                    .arg(this->nCoordGroups()), CODELOC );
    }

    QVarLengthArray<bool> to_remove(this->nCoordGroups());
    
    for (int i=0; i<this->nCoordGroups(); ++i)
    {
        to_remove[i] = false;
    }

    for (int i=int(start); i<int(start)+count; ++i)
    {
        to_remove[i] = true;
    }

    this->pvt_remove(to_remove);
}

/** Remove the ith CoordGroup from the array */
void CoordGroupArray::remove(quint32 i)
{
    this->remove(i, 1);
}

/** Remove the specified CoordGroups from the array */
void CoordGroupArray::remove(const QVarLengthArray<quint32> &idxs)
{
    if (idxs.isEmpty())
        return;

    QVarLengthArray<bool> to_remove(this->nCoordGroups());
    
    for (int i=0; i<this->nCoordGroups(); ++i)
    {
        to_remove[i] = false;
    }
    
    int nremoved = 0;
    
    for (int i=0; i<idxs.count(); ++i)
    {
        if (int(idxs[i]) >= this->nCoordGroups())
            throw SireError::invalid_index( QObject::tr(
                    "The list of indexes contains an invalid index (%1). "
                    "The number of CoordGroups is %2.")
                        .arg(idxs[i]).arg(this->nCoordGroups()), CODELOC );
                     
        if (not to_remove[i])
        {
            to_remove[i] = false;
            ++nremoved;
        }
    }
    
    if (nremoved == this->nCoordGroups())
    {
        this->operator=( CoordGroupArray() );
    }
    else
    {
        this->pvt_remove(to_remove);
    }
}

/** Update the CoordGroup at index i so that it has coordinates 'coords'
    (there are 'ncoords' coordinates in this array)
    
    There must contain the same number of coordinates as the existing
    CoordGroup at this index
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void CoordGroupArray::update(quint32 i, const Vector *coords, int ncoords)
{
    this->assertValidCoordGroup(i);
    
    CoordGroup &this_cgroup = d->cGroupData()[i];
    
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );
    
    if (cgdata->nCoords() != quint32(ncoords))
        throw SireError::incompatible_error( QObject::tr(
            "Cannot update the CoordGroup at index %1 as this CoordGroup "
            "has a different number of coordinates (%2) than are "
            "provided (%3).")
                .arg(i).arg(cgdata->nCoords()).arg(ncoords),
                    CODELOC );
                    
    //copy the coordinates
    Vector *old_coords = cgdata->coordsData();
    
    Vector *output = quickCopy<Vector>(old_coords, coords, ncoords);
                            
    BOOST_ASSERT(output == old_coords);
    
    //update the AABox
    *(cgdata->aaBox()) = AABox(coords, ncoords);
}

/** Update the CoordGroup at index i so that it has coordinates 'coords'
    
    There must contain the same number of coordinates as the existing
    CoordGroup at this index
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void CoordGroupArray::update(quint32 i, const QVector<Vector> &coords)
{
    this->update(i, coords.constData(), coords.count());
}

/** Update the CoordGroup at index i so that it is equal to 'cgroup'. Note
    that 'cgroup' must contain the same number of coordinates as the existing
    CoordGroup at this index
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void CoordGroupArray::update(quint32 i, const CoordGroup &cgroup)
{
    this->assertValidCoordGroup(i);
    
    CoordGroup &this_cgroup = d->cGroupData()[i];
    
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );
    
    if (cgdata->nCoords() != quint32(cgroup.count()))
        throw SireError::incompatible_error( QObject::tr(
            "Cannot update the CoordGroup at index %1 as this CoordGroup "
            "has a different number of coordinates (%2) than the "
            "new CoordGroup (%3).")
                .arg(i).arg(cgdata->nCoords()).arg(cgroup.count()),
                    CODELOC );
                    
    //copy the coordinates
    Vector *coords = cgdata->coordsData();
    
    Vector *output = quickCopy<Vector>(coords, cgroup.constData(),
                                       cgdata->nCoords());
                            
    BOOST_ASSERT(output == coords);
    
    //copy the AABox
    *(cgdata->aaBox()) = cgroup.aaBox();
}

/** Translate all of the coordinates in this array by 'delta' */
void CoordGroupArray::translate(const Vector &delta)
{
    if (delta.isZero())
        return;
        
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] += delta;
    }
    
    AABox *aaboxes = d->aaBoxData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i].translate(delta);
    }
}

/** Translate all of the coordinates in the ith CoordGroup by 'delta' */
void CoordGroupArray::translate(quint32 i, const Vector &delta)
{
    if (delta.isZero())
        return;
        
    this->assertValidCoordGroup(i);

    CoordGroup &this_cgroup = d->cGroupData()[i];
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] += delta;
    }
    
    cgdata->aaBox()->translate(delta);
}

/** Rotate all of coordinates in this array using the matrix 'rotmat'
    about the point 'point' */
void CoordGroupArray::rotate(const Matrix &rotmat, const Vector &point)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = SireMaths::rotate(coords[i], rotmat, point);
    }
    
    //now update all of the AABoxes
    CoordGroup *cgroups = d->cGroupData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        *( const_cast<CGData*>(cgroups[i].d.constData())->aaBox() ) 
               = AABox(cgroups[i]);
    }
}

/** Transform all of coordinates in this array using the transformation 't' */
void CoordGroupArray::transform(const Transform &t)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    coords = t.apply(coords, n);
    
    //now update all of the AABoxes
    CoordGroup *cgroups = d->cGroupData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        *( const_cast<CGData*>(cgroups[i].d.constData())->aaBox() ) 
               = AABox(cgroups[i]);
    }
}

/** Rotate all of the coordinates in this array using the quaternion
    'quat' around the point 'point' */
void CoordGroupArray::rotate(const Quaternion &quat, const Vector &point)
{
    this->rotate(quat.toMatrix(), point);
}

/** Rotate all of the coordinates in the CoordGroup at index 'i' using 
    the matrix 'rotmat' about the point 'point'
    
    \throw SireError::invalid_index
*/
void CoordGroupArray::rotate(quint32 i, const Matrix &rotmat, const Vector &point)
{
    this->assertValidCoordGroup(i);

    CoordGroup &this_cgroup = d->cGroupData()[i];
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = SireMaths::rotate(coords[i], rotmat, point);
    }
    
    *(cgdata->aaBox()) = AABox(this_cgroup);
}

/** Transform all of the coordinates in the CoordGroup at index 'i' using
    the transformation 't' 
 
    \throw SireError::invalid_index
*/
void CoordGroupArray::transform(quint32 i, const Transform &t)
{
    this->assertValidCoordGroup(i);

    CoordGroup &this_cgroup = d->cGroupData()[i];
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();
    
    coords = t.apply(coords, n);
    
    *(cgdata->aaBox()) = AABox(this_cgroup);
}

/** Rotate all of the coordinates in the CoordGroup at index 'i' using 
    the quaternion 'quat' about the point 'point'
    
    \throw SireError::invalid_index
*/
void CoordGroupArray::rotate(quint32 i, const Quaternion &quat, const Vector &point)
{
    this->rotate(i, quat.toMatrix(), point);
}

/** Map all of the coordinates in this array into the coordinate
    frame represented by 'axes' */
void CoordGroupArray::mapInto(const AxisSet &axes)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = axes.fromIdentity(coords[i]);
    }
    
    CoordGroup *cgroups = d->cGroupData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        *( const_cast<CGData*>(cgroups[i].d.constData())->aaBox() ) 
               = AABox(cgroups[i]);
    }
}

/** Map all of the coordinates of the CoordGroup at index 'i'
    into the coordinate frame represented by 'axes'
    
    \throw SireError::invalid_index
*/
void CoordGroupArray::mapInto(quint32 i, const AxisSet &axes)
{
    this->assertValidCoordGroup(i);

    CoordGroup &this_cgroup = d->cGroupData()[i];
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = axes.fromIdentity(coords[i]);
    }
    
    *(cgdata->aaBox()) = AABox(this_cgroup);
}

/** Change all of the coordinates in this array from the 
    coordinate frame 'from_frame' to the coordinate frame 'to_frame' */
void CoordGroupArray::changeFrame(const AxisSet &from_frame, 
                                  const AxisSet &to_frame)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = from_frame.toFrame(to_frame, coords[i]);
    }

    //update the AABoxes
    CoordGroup *cgroups = d->cGroupData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        *( const_cast<CGData*>(cgroups[i].d.constData())->aaBox() ) 
               = AABox(cgroups[i]);
    }

}

/** Change all of the coordinates in the ith CoordGroup from
    the coordinate frame 'from_frame' to the coordinate frame
    'to_frame'
    
    \throw SireError::invalid_index
*/
void CoordGroupArray::changeFrame(quint32 i, const AxisSet &from_frame, 
                                             const AxisSet &to_frame)
{
    this->assertValidCoordGroup(i);

    CoordGroup &this_cgroup = d->cGroupData()[i];
    CGData *cgdata = const_cast<CGData*>( this_cgroup.d.constData() );

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = from_frame.toFrame(to_frame, coords[i]);
    }
    
    *(cgdata->aaBox()) = AABox(this_cgroup);
}

const char* CoordGroupArray::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CoordGroupArray>() );
}

////////
//////// Implementation of CoordGroupArrayArray
////////

static const RegisterMetaType<CoordGroupArrayArray> r_cgarrayarray(NO_ROOT);

QDataStream& operator<<(QDataStream &ds,
                        const CGSharedPtr<CGArrayArrayData> &d)
{
    const CGArrayArrayData &array = *d;

    //serialise the prime metadata
    ds << quint32( array.nCGArrays() ) 
       << quint32( array.nCGroups() )
       << quint32( array.nCoords() );

    //serialise the number of CoordGroups in each CoordGroupArray
    for (quint32 i=0; i<array.nCGArrays(); ++i)
    {
        ds << quint32(array.cgArrayDataData()[i].nCGroups());
    }
    
    //serialise the number of coordinates in each CoordGroup
    for (quint32 i=0; i<array.nCGroups(); ++i)
    {
        ds << quint32(array.cgDataData()[i].nCoords());
    }
    
    //now the AABoxes
    for (quint32 i=0; i<array.nCGroups(); ++i)
    {
        ds << array.aaBoxData()[i];
    }
    
    //finally the coordinates
    for (quint32 i=0; i<array.nCoords(); ++i)
    {
        ds << array.coordsData()[i];
    }
    
    return ds;
}

static CGSharedPtr<CGArrayArrayData> createCGArrayArray(quint32 narrays,
                                        quint32 ngroups, quint32 ncoords)
{
    if (narrays == 0)
        return ::getSharedNull();

    //construct space for narrays arrays of ngroups CoordGroups of ncoords coordinates
    char *storage = CGMemory::create(narrays, ngroups, ncoords);
        
    CGArrayArrayData *array = (CGArrayArrayData*)storage;
    
    return CGSharedPtr<CGArrayArrayData>(array);
}

QDataStream& operator>>(QDataStream &ds,
                        CGSharedPtr<CGArrayArrayData> &d)
{
    quint32 ncgarrays, ncgroups, ncoords;
    
    //read the prime metadata
    ds >> ncgarrays >> ncgroups >> ncoords;

    //create space for this data
    CGSharedPtr<CGArrayArrayData> d2 = createCGArrayArray(ncgarrays, ncgroups,
                                                          ncoords);
        
    CGArrayArrayData &array = *d2;
                                                                                                                                            
    //dimension the arrays..
    for (quint32 i=0; i<ncgarrays; ++i)
    {
        quint32 n;
        ds >> n;
        
        array.setNCGroupsInArray(i, n);
    }

    for (quint32 i=0; i<ncgroups; ++i)
    {
        quint32 n;
        ds >> n;
        
        array.setNPointsInCGroup(i, n);
    }
    
    array.close();
    
    //copy in the AABoxes...
    for (quint32 i=0; i<ncgroups; ++i)
    {
        ds >> array.aaBoxData()[i];
    }
    
    //copy in the coordinates
    for (quint32 i=0; i<ncoords; ++i)
    {
        ds >> array.coordsData()[i];
    }
    
    d = d2;

    return ds;
}


/** This is a helper class that helps the loading and saving of a
    CGSharedPtr<CGArrayArrayData> */
struct GetCGSharedDataPointer
{
    static bool isEmpty(const CGSharedPtr<CGArrayArrayData> &d)
    {
        return d.constData() == 0;
    }

    static const void* value(const CGSharedPtr<CGArrayArrayData> &d)
    {
        return d.constData();
    }
    
    static void load(QDataStream &ds, CGSharedPtr<CGArrayArrayData> &d)
    {
        ds >> d;
    }
    
    static void save(QDataStream &ds, const CGSharedPtr<CGArrayArrayData> &d)
    {
        //null pointers will have been caught before here
        ds << d;
    }
};


SharedDataStream& operator<<(SharedDataStream &sds,
                             const CGSharedPtr<CGArrayArrayData> &d)
{
    sds.sharedSave< CGSharedPtr<CGArrayArrayData>, GetCGSharedDataPointer >(d);
    
    return sds;
}

SharedDataStream& operator>>(SharedDataStream &sds,
                             CGSharedPtr<CGArrayArrayData> &d)
{
    sds.sharedLoad< CGSharedPtr<CGArrayArrayData>, GetCGSharedDataPointer >(d);

    return sds;
}

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const CoordGroupArrayArray &cgarray)
{
    writeHeader(ds, r_cgarrayarray, 2);

    SharedDataStream sds(ds);
    
    sds << cgarray.d;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                       CoordGroupArrayArray &cgarray)
{
    VersionID v = readHeader(ds, r_cgarrayarray);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        sds >> cgarray.d;
    }
    else if (v == 1)
    {
        ds >> cgarray.d;
    }
    else
        throw version_error(v, "1,2", r_cgarrayarray, CODELOC);
    
    return ds;
}

/** Construct an empty array */
CoordGroupArrayArray::CoordGroupArrayArray()
                     : d( ::getSharedNull() )
{}

/** Construct an array that contains just the CoordGroup 'cgroup'
    (in a single array that contains this group) */
CoordGroupArrayArray::CoordGroupArrayArray(const CoordGroup &cgroup)
                     : d( cgroup.d->extract() )
{}

/** Construct an array that contains just the CoordGroupArray 'cgarray' */
CoordGroupArrayArray::CoordGroupArrayArray(const CoordGroupArray &cgarray)
                     : d( cgarray.d->extract() )
{}

/** Construct an array from a vector of arrays */
CoordGroupArrayArray::CoordGroupArrayArray(const QVector<CoordGroupArray> &cgarrays)
{
    //count the number of arrays, groups and coords
    quint32 narrays = cgarrays.count();
    quint32 ngroups = 0;
    quint32 ncoords = 0;
    
    const CoordGroupArray *cgroups_array = cgarrays.constData();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        const CoordGroupArray &cgarray = cgroups_array[i];
        
        quint32 this_ngroups = cgarray.count();
        ngroups += this_ngroups;
        
        const CoordGroup *cgarray_array = cgarray.constData();
        
        for (quint32 j=0; j<this_ngroups; ++j)
        {
            ncoords += cgarray_array[j].count();
        }
    }
    
    //create space for the contents...
    d = ::createCGArrayArray(narrays, ngroups, ncoords);
    
    CGArrayArrayData *dptr = d.data();
    
    quint32 igroup = 0;
    
    //now tell the space how the arrays are distributed
    for (quint32 i=0; i<narrays; ++i)
    {
        const CoordGroupArray &cgarray = cgroups_array[i];
        
        quint32 this_ngroups = cgarray.count();
        const CoordGroup *cgarray_array = cgarray.constData();
        
        dptr->setNCGroupsInArray(i, this_ngroups);
        
        for (quint32 j=0; j<this_ngroups; ++j)
        {
            dptr->setNPointsInCGroup(igroup, cgarray_array[j].count());
            ++igroup;
        }
    }
    
    dptr->close();
    
    //finally we must copy all of the data
    Vector *coords = dptr->coordsData();
    AABox *aaboxes = dptr->aaBoxData();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        //copy all of the coordinates in one block
        const CoordGroupArray &cgarray = cgroups_array[i];
        
        if (cgarray.nCoords() > 0)
        {
            Vector *output = quickCopy<Vector>( coords, cgarray.coordsData(),
                                                cgarray.nCoords() );
                                     
            BOOST_ASSERT(output == coords);
            
            coords += cgarray.nCoords();
        }
        
        //copy all of the AABoxes in one block
        if (cgarray.nCoordGroups() > 0)
        {
            AABox *output = quickCopy<AABox>( aaboxes, cgarray.aaBoxData(),
                                              cgarray.nCoordGroups() );
                                     
            BOOST_ASSERT(output == aaboxes);
            
            aaboxes += cgarray.nCoordGroups();
        }
    }
}

/** Construct from an array of array of CoordGroups */
CoordGroupArrayArray::CoordGroupArrayArray(
                        const QVector< QVector<CoordGroup> > &cgarrays)
{
    //count the number of arrays, groups and coords
    quint32 narrays = cgarrays.count();
    quint32 ngroups = 0;
    quint32 ncoords = 0;
    
    const QVector<CoordGroup> *cgroups_array = cgarrays.constData();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        const QVector<CoordGroup> &cgarray = cgroups_array[i];
        
        quint32 this_ngroups = cgarray.count();
        ngroups += this_ngroups;
        
        const CoordGroup *cgarray_array = cgarray.constData();
        
        for (quint32 j=0; j<this_ngroups; ++j)
        {
            ncoords += cgarray_array[j].count();
        }
    }
    
    //create space for the contents...
    d = ::createCGArrayArray(narrays, ngroups, ncoords);
    
    CGArrayArrayData *dptr = d.data();
    
    quint32 igroup = 0;
    
    //now tell the space how the arrays are distributed
    for (quint32 i=0; i<narrays; ++i)
    {
        const QVector<CoordGroup> &cgarray = cgroups_array[i];
        
        quint32 this_ngroups = cgarray.count();
        const CoordGroup *cgarray_array = cgarray.constData();
        
        dptr->setNCGroupsInArray(i, this_ngroups);
        
        for (quint32 j=0; j<this_ngroups; ++j)
        {
            dptr->setNPointsInCGroup(igroup, cgarray_array[j].count());
            ++igroup;
        }
    }
    
    dptr->close();
    
    //finally we must copy all of the data
    Vector *coords = dptr->coordsData();
    AABox *aabox = dptr->aaBoxData();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        const QVector<CoordGroup> &cgarray = cgroups_array[i];
        
        quint32 ngroups = cgarray.count();
        const CoordGroup *cgarray_array = cgarray.constData();
        
        for (quint32 j=0; j<ngroups; ++j)
        {
            const CoordGroup &cgroup = cgarray_array[j];
        
            if (cgroup.count() > 0)
            {
                Vector *output = quickCopy<Vector>( coords, cgroup.constData(),
                                                    cgroup.count() );
                                     
                BOOST_ASSERT(output == coords);
            
                coords += cgroup.count();
            }
        
            *aabox = cgroup.aaBox();
            ++aabox;
        }
    }
}

/** Construct from a triple-vector */
CoordGroupArrayArray::CoordGroupArrayArray(
                        const QVector< QVector< QVector<Vector> > > &points)
{
    if (points.isEmpty())
        return;

    QVector<CoordGroupArray> arrays(points.count());
    
    for (int i=0; i<points.count(); ++i)
    {
        arrays[i] = CoordGroupArray(points[i]);
    }
    
    this->operator=( CoordGroupArrayArray(arrays) );
}

/** Copy constructor */
CoordGroupArrayArray::CoordGroupArrayArray(const CoordGroupArrayArray &other)
                     : d(other.d)
{}

/** Destructor */
CoordGroupArrayArray::~CoordGroupArrayArray()
{}

/** Copy assignment operator */
CoordGroupArrayArray& CoordGroupArrayArray::operator=(
                                        const CoordGroupArrayArray &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool CoordGroupArrayArray::operator==(const CoordGroupArrayArray &other) const
{
    if (this->coordGroupData() == other.coordGroupData())
        return true;
        
    for (quint32 i=0; i<d->nCGroups(); ++i)
    {
        if (this->coordGroupData()[i] != other.coordGroupData()[i])
            return false;
    }
    
    return true;
}                                   

/** Comparison operator */
bool CoordGroupArrayArray::operator!=(const CoordGroupArrayArray &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this array of arrays */
QString CoordGroupArrayArray::toString() const
{
    if (this->count() == 0)
        return QObject::tr("CoordGroupArrayArray::null");
        
    QStringList lines;
    
    for (int i=0; i<this->count(); ++i)
    {
        lines.append( this->constData()[i].toString() );
    }
    
    return QObject::tr("CoordGroupArrayArray( aaBox() == %1\n%2\n)")
                .arg(this->aaBox().toString(), lines.join("\n"));
}

/** Merge this array of array of CoordGroups back into a single CoordGroup */
CoordGroup CoordGroupArrayArray::merge() const
{
    return CoordGroup(*this);
}

/** Assert that i is a valid index for a CoordGroupArray 

    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::assertValidCoordGroupArray(quint32 i) const
{
    if (i >= d->nCGArrays())
        throw SireError::invalid_index( QObject::tr(
            "Cannot access the CoordGroupArray at index %1. The "
            "number of CoordGroupArrays in this array is %2.")
                .arg(i).arg(d->nCGArrays()), CODELOC );
}

/** Assert that i is a valid index for a CoordGroupArray 

    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::assertValidIndex(quint32 i) const
{
    this->assertValidCoordGroupArray(i);
}

/** Assert that i is a valid index for a CoordGroup

    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::assertValidCoordGroup(quint32 i) const
{
    if (i >= d->nCGroups())
        throw SireError::invalid_index( QObject::tr(
            "Cannot access the CoordGroup at index %1. The "
            "number of CoordGroups in this array is %2.")
                .arg(i).arg(d->nCGroups()), CODELOC );
}

/** Assert that i is a valid index for a coordinate

    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::assertValidCoordinate(quint32 i) const
{
    if (i >= d->nCoords())
        throw SireError::invalid_index( QObject::tr(
            "Cannot access the coordinate at index %1. The "
            "number of coordinates in this array is %2.")
                .arg(i).arg(d->nCoords()), CODELOC );
}

/** Return the ith CoordGroupArray in this array

    \throw SireError::invalid_index
*/
const CoordGroupArray& CoordGroupArrayArray::operator[](quint32 i) const
{
    assertValidCoordGroupArray(i);
    return d->cgArrayData()[i];
}

/** Return the ith CoordGroupArray in this array

    \throw SireError::invalid_index
*/
const CoordGroupArray& CoordGroupArrayArray::at(quint32 i) const
{
    return this->operator[](i);
}

/** Return the number of CoordGroupArrays in this array */
int CoordGroupArrayArray::count() const
{
    return d->nCGArrays();
}

/** Return the number of CoordGroupArrays in this array */
int CoordGroupArrayArray::size() const
{
    return d->nCGArrays();
}

/** Return the number of CoordGroupArrays in this array */
int CoordGroupArrayArray::nCoordGroupArrays() const
{
    return d->nCGArrays();
}

/** Return the total number of CoordGroups in all of the CoordGroupArrays
    in this array */
int CoordGroupArrayArray::nCoordGroups() const
{
    return d->nCGroups();
}

/** Return the total number of points in all of the CoordGroups in
    all of the CoordGroupArrays in this array */
int CoordGroupArrayArray::nCoords() const
{
    return d->nCoords();
}

/** Return an AABox that complete encompasses all of the CoordGroups 
    in this array */
AABox CoordGroupArrayArray::aaBox() const
{
    int ncgroups = this->nCoordGroups();

    if (ncgroups == 1)
        return this->coordGroupData()[0].aaBox();
    else if (ncgroups == 0)
        return AABox();
    else
    {
        const CoordGroup *coords_array = this->coordGroupData();
        
        Vector mincoords = coords_array[0].aaBox().minCoords();
        Vector maxcoords = coords_array[0].aaBox().maxCoords();
        
        for (int i=1; i<ncgroups; ++i)
        {
            mincoords.setMin(coords_array[i].aaBox().minCoords());
            maxcoords.setMax(coords_array[i].aaBox().maxCoords());
        }
        
        return AABox::from(mincoords, maxcoords);
    }
}

/** Return a raw pointer to the array of CoordGroupArrays */
const CoordGroupArray* CoordGroupArrayArray::data() const
{
    return d->cgArrayData();
}

/** Return a raw pointer to the array of CoordGroupArrays */
const CoordGroupArray* CoordGroupArrayArray::constData() const
{
    return d->cgArrayData();
}

/** Return a raw pointer to the array of CoordGroups */
const CoordGroup* CoordGroupArrayArray::coordGroupData() const
{
    return d->cGroupData();
}

/** Return a raw pointer to the array of CoordGroups */
const CoordGroup* CoordGroupArrayArray::constCoordGroupData() const
{
    return d->cGroupData();
}

/** Return a raw pointer to the array of coordinates */
const Vector* CoordGroupArrayArray::coordsData() const
{
    return d->coordsData();
}

/** Return a raw pointer to the array of coordinates */
const Vector* CoordGroupArrayArray::constCoordsData() const
{
    return d->coordsData();
}

/** Return a raw pointer to the array of AABoxes */
const AABox* CoordGroupArrayArray::aaBoxData() const
{
    return d->aaBoxData();
}

/** Return a raw pointer to the array of AABoxes */
const AABox* CoordGroupArrayArray::constAABoxData() const
{
    return d->aaBoxData();
}

/** Update the ith CoordGroupArray so that its contents are equal to 'array'. 
    Note that you can only do this if the configuration of the
    CoordGroupArray at index i is identical to 'array' (same number
    of CoordGroups, with the same number of coordinates in 
    equivalent CoordGroups)
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void CoordGroupArrayArray::update(quint32 i, const CoordGroupArray &array)
{
    this->assertValidCoordGroupArray(i);

    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    quint32 n = cgarraydata->nCoords();

    if (n != quint32(array.count()))
        throw SireError::incompatible_error( QObject::tr(
            "Cannot update the CoordGroupArray at index %1 as "
            "it contains %2 CoordGroups, but the new array has %3.")
                .arg(i).arg(n).arg(array.count()), CODELOC );

    //ensure all of the CoordGroups have the same number of points...
    const CGArrayData *new_cgarray = array.d.constData();
    
    for (quint32 j=0; j<n; ++j)
    {
        if (cgarraydata->cGroupData()[j].count() 
                != new_cgarray->cGroupData()[j].count())
        {
            throw SireError::incompatible_error( QObject::tr(
                "Cannot update the CoordGroupArray at index %1 as "
                "the number of coordinates in the CoordGroup at index "
                "%2 (%3) is not the same as the number of "
                "CoordGroups at the same index in 'array' (%4).")
                    .arg(i).arg(j)
                    .arg(cgarraydata->cGroupData()[j].count())
                    .arg(new_cgarray->cGroupData()[j].count()), CODELOC );
        }
    }
    
    //ok, we can now copy the AABoxes as a block, then the 
    //CoordGroups as a block
    AABox *output = quickCopy<AABox>( cgarraydata->aaBoxData(),
                                      new_cgarray->aaBoxData(),
                                      cgarraydata->nCGroups() );
                             
    BOOST_ASSERT( output == cgarraydata->aaBoxData() );
    
    Vector *out2 = quickCopy<Vector>( cgarraydata->coordsData(),
                                      new_cgarray->coordsData(),
                                      cgarraydata->nCoords() );
                       
    BOOST_ASSERT( out2 == cgarraydata->coordsData() );
}

/** Update the jth CoordGroup of the ith CoordGroupArray so that its
    coordinates are equal to 'cgroup'. Note that you can only do this if
    the number of coordinates in this CoordGroup is the same as
    in 'cgroup'
    
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
*/
void CoordGroupArrayArray::update(quint32 i, quint32 j, const CoordGroup &cgroup)
{
    this->assertValidCoordGroupArray(i);

    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    if (j >= cgarraydata->nCGroups())
        throw SireError::invalid_index( QObject::tr(
            "There is no CoordGroup at index %1 in the array at "
            "index %2. The number of CoordGroups in this array "
            "is %3.")
                .arg(j).arg(i).arg(cgarraydata->nCGroups()), CODELOC );
    
    CGData *cgdata = const_cast<CGData*>( cgarraydata->cGroupData()[j]
                                            .d.constData() );
    
    if (cgdata->nCoords() != quint32(cgroup.count()))
        throw SireError::incompatible_error( QObject::tr(
            "Cannot update the CoordGroup at index %1 in the array "
            "at index %2, as the number of coordinates in this "
            "CoordGroup (%3) is not equal to the new CoordGroup (%4).")
                .arg(j).arg(i).arg(cgdata->nCoords()).arg(cgroup.count()),
                    CODELOC );

    //copy the coordinates and AABox in one go
    *(cgdata->aaBox()) = cgroup.aaBox();
    
    Vector *output = quickCopy<Vector>( cgdata->coordsData(),
                                        cgroup.constData(),
                                        cgdata->nCoords() );
                             
    BOOST_ASSERT( output == cgdata->coordsData() );
}

/** Translate all of the coordinates held in this container by 'delta' */
void CoordGroupArrayArray::translate(const Vector &delta)
{
    if (delta.isZero())
        return;
        
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] += delta;
    }
    
    //update all of the AABoxes
    AABox *aaboxes = d->aaBoxData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i].translate(delta);
    }
}

/** Translate all of the points in the ith array by delta

    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::translate(quint32 i, const Vector &delta)
{
    this->assertValidCoordGroupArray(i);

    if (delta.isZero())
        return;
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    quint32 n = cgarraydata->nCoords();
    Vector *coords = cgarraydata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] += delta;
    }
    
    n = cgarraydata->nCGroups();
    AABox *aaboxes = cgarraydata->aaBoxData();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i].translate(delta);
    }
}

/** Assert that there is a jth CoordGroup in the ith array 

    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::assertValidCoordGroup(quint32 i, quint32 j) const
{
    this->assertValidCoordGroupArray(i);
    
    if (j >= quint32(d->cgArrayData()[i].nCoordGroups()))
        throw SireError::invalid_index( QObject::tr(
            "There is no CoordGroup at index %1 in the CoordGroupArray "
            "at index %2. The number of CoordGroups in this array is %3.")
                .arg(j).arg(i).arg(d->cgArrayData()[i].nCoordGroups()), CODELOC );
}

/** Translate all of the points in the jth CoordGroup of the ith
    array by 'delta'
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::translate(quint32 i, quint32 j, const Vector &delta)
{
    this->assertValidCoordGroup(i,j);
    
    if (delta.isZero())
        return;
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    CGData *cgdata = const_cast<CGData*>(cgarraydata->cGroupData()[j].d.constData());

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] += delta;
    }
    
    cgdata->aaBox()->translate(delta);
}

/** Rotate all of the points in this container using the matrix 'rotmat'
    about the point 'point' */
void CoordGroupArrayArray::rotate(const Matrix &rotmat, const Vector &point)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = SireMaths::rotate(coords[i], rotmat, point);
    }
    
    //update all of the AABoxes
    AABox *aaboxes = d->aaBoxData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox(d.constData()->aaBoxData()[i]);
    }
}

/** Rotate all of the points in this container using the matrix 'rotmat'
    about the point 'point' */
void CoordGroupArrayArray::transform(const Transform &t)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    coords = t.apply(coords, n);
    
    //update all of the AABoxes
    AABox *aaboxes = d->aaBoxData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox(d.constData()->aaBoxData()[i]);
    }
}

/** Rotate all of the points in this container using the quaternion 'quat'
    about the point 'point' */
void CoordGroupArrayArray::rotate(const Quaternion &quat, const Vector &point)
{
    this->rotate(quat.toMatrix(), point);
}

/** Rotate all of the points in the ith array using the matrix 'rotmat'
    about the point 'point' 
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::rotate(quint32 i, const Matrix &rotmat, 
                                  const Vector &point)
{
    this->assertValidCoordGroupArray(i);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    quint32 n = cgarraydata->nCoords();
    Vector *coords = cgarraydata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = SireMaths::rotate(coords[i], rotmat, point);
    }
    
    n = cgarraydata->nCGroups();
    AABox *aaboxes = cgarraydata->aaBoxData();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox( cgarraydata->cGroupData()[i] );
    }
}

/** Transform all of the points in the ith array using the 
    transformation 't'
 
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::transform(quint32 i, const Transform &t)
{
    this->assertValidCoordGroupArray(i);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    quint32 n = cgarraydata->nCoords();
    Vector *coords = cgarraydata->coordsData();

    coords = t.apply(coords, n);
    
    n = cgarraydata->nCGroups();
    AABox *aaboxes = cgarraydata->aaBoxData();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox( cgarraydata->cGroupData()[i] );
    }
}

/** Rotate all of the points in the ith array using the quaternion 'quat'
    about the point 'point' 
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::rotate(quint32 i, const Quaternion &quat, 
                                  const Vector &point)
{
    this->rotate(i, quat.toMatrix(), point);
}

/** Rotate all of the points in the jth CoordGroup in the 
    ith array using the matrix 'rotmat' about the point 'point' 
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::rotate(quint32 i, quint32 j,
                                  const Matrix &rotmat, const Vector &point)
{
    this->assertValidCoordGroup(i,j);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    CoordGroup &cgroup = cgarraydata->cGroupData()[j];
    CGData *cgdata = const_cast<CGData*>(cgroup.d.constData());

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = SireMaths::rotate(coords[i], rotmat, point);
    }
    
    *(cgdata->aaBox()) = AABox(cgroup);
}

/** Transform all of the points in the jth CoordGroup in the
    ith array using the transformation 't'
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::transform(quint32 i, quint32 j, const Transform &t)
{
    this->assertValidCoordGroup(i,j);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    CoordGroup &cgroup = cgarraydata->cGroupData()[j];
    CGData *cgdata = const_cast<CGData*>(cgroup.d.constData());

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();

    coords = t.apply(coords, n);
    
    *(cgdata->aaBox()) = AABox(cgroup);
}

/** Rotate all of the points in the jth CoordGroup in the 
    ith array using the quaternion 'quat' about the point 'point' 
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::rotate(quint32 i, quint32 j,
                                  const Quaternion &quat, const Vector &point)
{
    this->rotate(i, j, quat.toMatrix(), point);
}

/** Map all of the points in this array into the coordinate frame
    represented by 'axes' */
void CoordGroupArrayArray::mapInto(const AxisSet &axes)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = axes.fromIdentity(coords[i]);
    }
    
    //update all of the AABoxes
    AABox *aaboxes = d->aaBoxData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox(d.constData()->aaBoxData()[i]);
    }
}

/** Map all of the points in the ith array into the coordinate frame
    represented by 'axes'
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::mapInto(quint32 i, const AxisSet &axes)
{
    this->assertValidCoordGroupArray(i);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    quint32 n = cgarraydata->nCoords();
    Vector *coords = cgarraydata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = axes.fromIdentity(coords[i]);
    }
    
    n = cgarraydata->nCGroups();
    AABox *aaboxes = cgarraydata->aaBoxData();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox( cgarraydata->cGroupData()[i] );
    }
}

/** Map all the points in the CoordGroup at index j in the array at 
    index i into the coordinate frame represented by 'axes'
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::mapInto(quint32 i, quint32 j, const AxisSet &axes)
{
    this->assertValidCoordGroup(i,j);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    CoordGroup &cgroup = cgarraydata->cGroupData()[j];
    CGData *cgdata = const_cast<CGData*>(cgroup.d.constData());

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = axes.fromIdentity(coords[i]);
    }
    
    *(cgdata->aaBox()) = AABox(cgroup);
}

/** Map all of the points in this array from the coordinate frame
    'from_frame' to the coordinate frame 'to_frame' */
void CoordGroupArrayArray::changeFrame(const AxisSet &from_frame, 
                                       const AxisSet &to_frame)
{
    Vector *coords = d->coordsData();
    quint32 n = d.constData()->nCoords();
    
    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = from_frame.toFrame(to_frame, coords[i]);
    }
    
    //update all of the AABoxes
    AABox *aaboxes = d->aaBoxData();
    n = d.constData()->nCGroups();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox(d.constData()->aaBoxData()[i]);
    }
}

/** Map all of the points in this array from the coordinate frame
    'from_frame' to the coordinate frame 'to_frame'
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::changeFrame(quint32 i, const AxisSet &from_frame, 
                                                  const AxisSet &to_frame)
{
    this->assertValidCoordGroupArray(i);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    quint32 n = cgarraydata->nCoords();
    Vector *coords = cgarraydata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = from_frame.toFrame(to_frame, coords[i]);
    }
    
    n = cgarraydata->nCGroups();
    AABox *aaboxes = cgarraydata->aaBoxData();
    
    for (quint32 i=0; i<n; ++i)
    {
        aaboxes[i] = AABox( cgarraydata->cGroupData()[i] );
    }
}

/** Map all of the points in the CoordGroup at index j in the
    array at index i from the coordinate frame 'from_frame' to 
    the coordinate frame 'to_frame'
    
    \throw SireError::invalid_index
*/
void CoordGroupArrayArray::changeFrame(quint32 i, quint32 j,
                                       const AxisSet &from_frame, 
                                       const AxisSet &to_frame)
{
    this->assertValidCoordGroup(i,j);
    
    CoordGroupArray &this_cgarray = d->cgArrayData()[i];
    CGArrayData *cgarraydata = const_cast<CGArrayData*>( this_cgarray.d.constData() );

    CoordGroup &cgroup = cgarraydata->cGroupData()[j];
    CGData *cgdata = const_cast<CGData*>(cgroup.d.constData());

    quint32 n = cgdata->nCoords();
    Vector *coords = cgdata->coordsData();

    for (quint32 i=0; i<n; ++i)
    {
        coords[i] = from_frame.toFrame(to_frame, coords[i]);
    }
    
    *(cgdata->aaBox()) = AABox(cgroup);
}

const char* CoordGroupArrayArray::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CoordGroupArrayArray>() );
}
