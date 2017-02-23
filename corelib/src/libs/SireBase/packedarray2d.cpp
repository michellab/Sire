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

#include "packedarray2d.h"
#include "packedarray2d.hpp"

#include "SireError/errors.h"
#include "SireStream/datastream.h"

using namespace SireBase;
using namespace SireBase::detail;

using namespace SireStream;

void SIREBASE_EXPORT 
SireBase::detail::throwCannotConvertVariantError(const char *this_type, 
                                                 const char *type_t,
                                                 const QString &codeloc)
{
    throw SireError::invalid_cast( QObject::tr(
            "Cannot cast a QVariant of type %1 so that is can be "
            "placed into a PackedArray2D of type %2.")
                .arg(this_type).arg(type_t), CODELOC );
}

void SIREBASE_EXPORT SireBase::detail::throwPackedArray2D_invalidIndex(
                                                    quint32 i, quint32 nvals)
{
    throw SireError::invalid_index( QObject::tr(
        "The index %1 is invalid, as the number of arrays equals %2.")
            .arg(i).arg(nvals), CODELOC );
}

void SIREBASE_EXPORT SireBase::detail::throwPackedArray2D_Array_invalidIndex(
                                                          quint32 i, quint32 nvals)
{
    throw SireError::invalid_index( QObject::tr(
        "The index %1 is invalid, as the number of objects in the array equals %2.")
            .arg(i).arg(nvals), CODELOC );
}

void SIREBASE_EXPORT SireBase::detail::throwPackedArray2D_Array_incompatibleError(
                                                    quint32 this_sz, quint32 other_sz)
{
    throw SireError::incompatible_error( QObject::tr(
        "Cannot update as the PackedArray<T>::Arrays are incompatible! "
        "(number of objects in this array is %1, while the number of "
        "objects in the other array is %2.")
            .arg(this_sz).arg(other_sz), CODELOC );
}       

static const RegisterMetaType<PackedArray2DDataBase> r_parray( MAGIC_ONLY, NO_ROOT,
                                                      "SireBase::PackedArray2D<T>" );

void SIREBASE_EXPORT 
SireBase::detail::writePackedArray2DHeader(QDataStream &ds)
{
    writeHeader(ds, r_parray, 2);
}

quint32 SIREBASE_EXPORT 
SireBase::detail::readPackedArray2DHeader(QDataStream &ds)
{
    VersionID v = readHeader(ds, r_parray);
    
    if (v != 1 and v != 2)
        throw version_error(v, "1,2", r_parray, CODELOC);

    return v;
}

static const RegisterMetaType<PackedArray2D_ArrayDataBase> r_parrayarray( MAGIC_ONLY, NO_ROOT,
                                                "SireBase::PackedArray2D<T>::Array" );

void SIREBASE_EXPORT 
SireBase::detail::writePackedArray2DArrayHeader(QDataStream &ds, quint32 version)
{
    writeHeader(ds, r_parrayarray, version);
}

void SIREBASE_EXPORT 
SireBase::detail::readPackedArray2DArrayHeader(QDataStream &ds, quint32 version)
{
    VersionID v = readHeader(ds, r_parrayarray);
    
    if (v != version)
        throw version_error(v, QString::number(version), r_parrayarray, CODELOC);
}

////////
//////// Implementation of PackedArray2DMemoryBase
////////

/** Return the number of bytes necessary to store the array of arrays */
quint32 PackedArray2DMemoryBase::getSize(quint32 narrays, quint32 nvalues, 
                                         quint32 sizeof_PackedArray2DData,
                                         quint32 sizeof_PackedArray2D_Array,
                                         quint32 sizeof_PackedArray2D_ArrayData,
                                         quint32 sizeof_T)
{
    //need to add space for null arrays
    if (narrays == 0)
        narrays = 1;

    return sizeof_PackedArray2DData +
           narrays * sizeof_PackedArray2D_Array +
           narrays * sizeof_PackedArray2D_ArrayData +
           nvalues * sizeof_T;
}

/** Allocate the memory necessary to hold the array of arrays */
char* PackedArray2DMemoryBase::create(quint32 narrays, quint32 nvalues, 
                                      quint32 sizeof_PackedArray2DData,
                                      quint32 sizeof_PackedArray2D_Array,
                                      quint32 sizeof_PackedArray2D_ArrayData,
                                      quint32 sizeof_T)
{
    quint32 sz = PackedArray2DMemoryBase::getSize(narrays, nvalues,
                                                  sizeof_PackedArray2DData,
                                                  sizeof_PackedArray2D_Array,
                                                  sizeof_PackedArray2D_ArrayData,
                                                  sizeof_T);
                                          
    return new char[sz];
}

/////////
///////// Implementation of PackedArray2DDataBase
/////////

/** Null constructor */
PackedArray2DDataBase::PackedArray2DDataBase()
                      : RefCountData(), array0(0), arraydata0(0),
                        narrays(0), value0(0), nvalues(0)
{}

/** Construct to hold the specified number of arrays and values */
PackedArray2DDataBase::PackedArray2DDataBase(quint32 na, quint32 nv)
                      : RefCountData(), array0(0), arraydata0(0),
                        narrays(na), value0(0), nvalues(nv)
{}

/** Copy constructor */
PackedArray2DDataBase::PackedArray2DDataBase(const PackedArray2DDataBase &other)
                      : RefCountData(), array0(other.array0),
                        arraydata0(other.arraydata0), narrays(other.narrays),
                      value0(other.value0), nvalues(other.nvalues)
{}

/** Destructor */
PackedArray2DDataBase::~PackedArray2DDataBase()
{}

/** Increse the reference count of this object */
void PackedArray2DDataBase::incref()
{
    this->ref.ref();
}

/** Return the number of arrays in this container */
quint32 PackedArray2DDataBase::nArrays() const
{
    return narrays;
}

/** Return the total number of objects in this container */
quint32 PackedArray2DDataBase::nValues() const
{
    return nvalues;
}

/** Return a raw pointer to the memory used to store this container */
char* PackedArray2DDataBase::memory()
{
    return PackedArray2DMemoryBase::getRoot( (char*)this, 0 );
}

/** Return a raw pointer to the memory used to store this container */
const char* PackedArray2DDataBase::memory() const
{
    return PackedArray2DMemoryBase::getRoot( (const char*)this, 0 );
}

//////////
////////// Implementation of PackedArray2D_ArrayDataBase
//////////

/** Null constructor */
PackedArray2D_ArrayDataBase::PackedArray2D_ArrayDataBase()
                            : this_array(0), value0(0), nvalues(0)
{}

/** Construct, telling this array where it is in the storage array */
PackedArray2D_ArrayDataBase::PackedArray2D_ArrayDataBase(quint32 this_idx)
                            : this_array(this_idx), value0(0), nvalues(0)
{}

/** Copy constructor */
PackedArray2D_ArrayDataBase::PackedArray2D_ArrayDataBase(
                                    const PackedArray2D_ArrayDataBase &other)
                            : this_array(other.this_array),
                              value0(other.value0), nvalues(other.nvalues)
{}

/** Destructor */
PackedArray2D_ArrayDataBase::~PackedArray2D_ArrayDataBase()
{}

/** Return the number of values in the array */
quint32 PackedArray2D_ArrayDataBase::nValues() const
{
    return nvalues;
}

/** Return a const pointer to the start of the storage
    array that contains this ArrayData */
const char* PackedArray2D_ArrayDataBase::memory() const
{
    return PackedArray2DMemoryBase::getRoot( (const char*)this, this_array );
}

/** Return a pointer to the start of the storage array that
    contains this ArrayData */
char* PackedArray2D_ArrayDataBase::memory()
{
    return PackedArray2DMemoryBase::getRoot( (char*)this, this_array );
}

/// fully instantiate key PackedArray2D classes
namespace SireBase
{
    template class PackedArray2D<QString>;
    template class PackedArray2D<qint64>;
    template class PackedArray2D<double>;
    template class PackedArray2D<QVariant>;
}
