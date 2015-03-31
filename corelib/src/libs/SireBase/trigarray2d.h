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

#ifndef SIREBASE_TRIGARRAY2D_H
#define SIREBASE_TRIGARRAY2D_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class TrigArray2DBase;
}

QDataStream& operator<<(QDataStream&, const SireBase::TrigArray2DBase&);
QDataStream& operator>>(QDataStream&, SireBase::TrigArray2DBase&);

namespace SireBase
{

/** Base class of the TrigArray2D<T> class

    @author Christopher Woods
*/
class SIREMOL_EXPORT TrigArray2DBase
{

friend QDataStream& ::operator<<(QDataStream&, const TrigArray2DBase&);
friend QDataStream& ::operator>>(QDataStream&, TrigArray2DBase&);

public:
    ~TrigArray2DBase();
    
    int size() const;
    int count() const;
    
    int nRows() const;
    int nColumns() const;
    
    int map(int i, int j) const;
    
    int offset(int i, int j) const;
    int checkedOffset(int i, int j) const;
    
    void assertValidIndex(int i, int j) const;

protected:
    TrigArray2DBase();

    TrigArray2DBase(int dimension);

    TrigArray2DBase(const TrigArray2DBase &other);

    TrigArray2DBase& operator=(const TrigArray2DBase &other);

    bool operator==(const TrigArray2DBase &other) const;
    bool operator!=(const TrigArray2DBase &other) const;

private:
    void throwInvalidIndex(int i, int j) const;
    
    /** The dimension of this square matrix */
    qint32 dim;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the number of rows in the matrix */
inline int TrigArray2DBase::nRows() const
{
    return dim;
}

/** Return the number of columns in the matrix */
inline int TrigArray2DBase::nColumns() const
{
    return dim;
}

/** Return the number of unique items in this array (the size of 
    the underlying 1D array) */
inline int TrigArray2DBase::size() const
{
	return (dim*dim + dim)/2;
}

/** Return the number of unique items in this array (the size of 
    the underlying 1D array) */
inline int TrigArray2DBase::count() const
{
	return TrigArray2DBase::size();
}

/** Return the offset into the array of the value at index [i,j] */
inline int TrigArray2DBase::offset(int i, int j) const
{
    if (i <= j)
        return (2*(j + i*dim) - i - i*i) / 2;
    else
        return (2*(i + j*dim) - j - j*j) / 2;
}

/** Map the 2D index (i,j) into the 1D index into memory */
inline int TrigArray2DBase::map(int i, int j) const
{
    return TrigArray2DBase::offset(i, j);
}

/** Assert that the index (i,j) is valid for this matrix

    \throw SireError::invalid_index
*/
inline void TrigArray2DBase::assertValidIndex(int i, int j) const
{
    if (i < 0 or i >= dim or j < 0 or j >= dim)
        throwInvalidIndex(i,j);
}

inline int TrigArray2DBase::checkedOffset(int i, int j) const
{
    TrigArray2DBase::assertValidIndex(i,j);
    return TrigArray2DBase::offset(i,j);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireBase::TrigArray2DBase )

SIRE_END_HEADER

#endif
