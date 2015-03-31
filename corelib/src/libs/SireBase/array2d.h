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

#ifndef SIREBASE_ARRAY2D_H
#define SIREBASE_ARRAY2D_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class Array2DBase;
}

QDataStream& operator<<(QDataStream&, const SireBase::Array2DBase&);
QDataStream& operator>>(QDataStream&, SireBase::Array2DBase&);

namespace SireBase
{

/** Base class of the Array2D<T> class

    @author Christopher Woods
*/
class SIREMOL_EXPORT Array2DBase
{

friend QDataStream& ::operator<<(QDataStream&, const Array2DBase&);
friend QDataStream& ::operator>>(QDataStream&, Array2DBase&);

public:
    ~Array2DBase();
    
    int nRows() const;
    int nColumns() const;
    
    int offset(int i, int j) const;
    int checkedOffset(int i, int j) const;
    
    int map(int i, int j) const;
    
    void assertValidIndex(int i, int j) const;

protected:
    Array2DBase();

    Array2DBase(int nrows, int ncolumns);

    Array2DBase(const Array2DBase &other);

    Array2DBase& operator=(const Array2DBase &other);

    bool operator==(const Array2DBase &other) const;
    bool operator!=(const Array2DBase &other) const;

private:
    void throwInvalidIndex(int i, int j) const;
    
    /** The number of rows and columns in the matrix */
    qint32 nrows, ncolumns;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the number of rows in the matrix */
inline int Array2DBase::nRows() const
{
    return nrows;
}

/** Return the number of columns in the matrix */
inline int Array2DBase::nColumns() const
{
    return ncolumns;
}

/** Return the location in the 1D array of the item at index [i,j] */
inline int Array2DBase::offset(int i, int j) const
{
    return (i*ncolumns) + j;
}

/** Map the 2D index (i,j) into the 1D index into memory */
inline int Array2DBase::map(int i, int j) const
{
    return Array2DBase::offset(i,j);
}

/** Assert that the index (i,j) is valid for this matrix

    \throw SireError::invalid_index
*/
inline void Array2DBase::assertValidIndex(int i, int j) const
{
    if (i < 0 or i >= nrows or j < 0 or j >= ncolumns)
        throwInvalidIndex(i,j);
}

/** Return the location in the 1D array of the item at index [i,j] 

    \throw SireError::invalid_index
*/
inline int Array2DBase::checkedOffset(int i, int j) const
{
    Array2DBase::assertValidIndex(i,j);
    return Array2DBase::offset(i,j);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireBase::Array2DBase )

SIRE_END_HEADER

#endif
