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

#include "trigarray2d.h"
#include "trigarray2d.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<TrigArray2DBase> r_array2d( MAGIC_ONLY, NO_ROOT,
                                                          "SireBase::TrigArray2D<T>" );
                                                      
/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds,
                                        const TrigArray2DBase &array2d)
{
    writeHeader(ds, r_array2d, 1);
    ds << array2d.dim;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds,
                                        TrigArray2DBase &array2d)
{
    VersionID v = readHeader(ds, r_array2d);
    
    if (v == 1)
    {
        ds >> array2d.dim;
    }
    else
        throw version_error(v, "1", r_array2d, CODELOC);
    
    return ds;
}

/** Null constructor */
TrigArray2DBase::TrigArray2DBase() : dim(0)
{}

/** Construct with the specified dimension [dim,dim] */
TrigArray2DBase::TrigArray2DBase(int dimension)
            : dim(dimension)
{
    if (dimension < 0)
        dim = 0;
}

/** Copy constructor */
TrigArray2DBase::TrigArray2DBase(const TrigArray2DBase &other)
                : dim(other.dim)
{}

/** Destructor */
TrigArray2DBase::~TrigArray2DBase()
{}

/** Copy assignment operator */
TrigArray2DBase& TrigArray2DBase::operator=(const TrigArray2DBase &other)
{
    dim = other.dim;
    return *this;
}

/** Comparison operator */
bool TrigArray2DBase::operator==(const TrigArray2DBase &other) const
{
    return dim == other.dim;
}

/** Comparison operator */
bool TrigArray2DBase::operator!=(const TrigArray2DBase &other) const
{
    return dim != other.dim;
}

/** Throw an invalid index exception */
void TrigArray2DBase::throwInvalidIndex(int i, int j) const
{
    throw SireError::invalid_index( QObject::tr(
        "Index (%1,%2) is not valid for this matrix. The "
        "number of rows is %3, and number of columns is %4.")
            .arg(i).arg(j).arg(dim).arg(dim), CODELOC );
}
