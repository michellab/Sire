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

#include "array2d.h"
#include "array2d.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Array2DBase> r_array2d( MAGIC_ONLY, NO_ROOT,
                                                      "SireBase::Array2D<T>" );
                                                      
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const Array2DBase &array2d)
{
    writeHeader(ds, r_array2d, 2);
    ds << array2d.nrows << array2d.ncolumns;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        Array2DBase &array2d)
{
    VersionID v = readHeader(ds, r_array2d);
    
    if (v == 2)
    {
        ds >> array2d.nrows >> array2d.ncolumns;
    }
    else if (v == 1)
    {
        quint32 nrows, ncolumns;
        ds >> nrows >> ncolumns;
        
        array2d.nrows = nrows;
        array2d.ncolumns = ncolumns;
    }
    else
        throw version_error(v, "1,2", r_array2d, CODELOC);
    
    return ds;
}

/** Null constructor */
Array2DBase::Array2DBase() : nrows(0), ncolumns(0)
{}

/** Construct with the specified number of rows and columns */
Array2DBase::Array2DBase(int nr, int nc)
            : nrows(nr), ncolumns(nc)
{
    if (nr < 0)
        nrows = 0;
        
    if (nc < 0)
        ncolumns = 0;
}

/** Copy constructor */
Array2DBase::Array2DBase(const Array2DBase &other)
            : nrows(other.nrows), ncolumns(other.ncolumns)
{}

/** Destructor */
Array2DBase::~Array2DBase()
{}

/** Copy assignment operator */
Array2DBase& Array2DBase::operator=(const Array2DBase &other)
{
    nrows = other.nrows;
    ncolumns = other.ncolumns;
    return *this;
}

/** Comparison operator */
bool Array2DBase::operator==(const Array2DBase &other) const
{
    return nrows == other.nrows and ncolumns == other.ncolumns;
}

/** Comparison operator */
bool Array2DBase::operator!=(const Array2DBase &other) const
{
    return nrows != other.nrows or ncolumns != other.ncolumns;
}

/** Throw an invalid index exception */
void Array2DBase::throwInvalidIndex(int i, int j) const
{
    throw SireError::invalid_index( QObject::tr(
        "Index (%1,%2) is not valid for this matrix. The "
        "number of rows is %3, and number of columns is %4.")
            .arg(i).arg(j).arg(nrows).arg(ncolumns), CODELOC );
}
