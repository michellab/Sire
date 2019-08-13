/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREVOL_GRIDINFO_H
#define SIREVOL_GRIDINFO_H

#include "SireVol/aabox.h"
#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class GridIndex;
class GridInfo;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::GridIndex&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::GridIndex&);

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::GridInfo&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::GridInfo&);

namespace SireMaths
{
class MultiFloat;
class MultiInt;
}

namespace SireVol
{

using SireMaths::Vector;
using SireVol::AABox;
using SireUnits::Dimension::Length;

using SireMaths::MultiInt;
using SireMaths::MultiFloat;

/** Very simple class providing a grid index */
class SIREVOL_EXPORT GridIndex
{
public:
    GridIndex(int i=0, int j=0, int k=0) : _i(i), _j(j), _k(k)
    {
        if (i < 0 or j < 0 or k < 0)
        {
            _i = 0;
            _j = 0;
            _k = -1;
        }
    }
    
    GridIndex(const GridIndex &other) : _i(other._i), _j(other._j), _k(other._k)
    {}
    
    ~GridIndex()
    {}
    
    GridIndex& operator=(const GridIndex &other)
    {
        _i = other._i;
        _j = other._j;
        _k = other._k;
        return *this;
    }
    
    bool operator==(const GridIndex &other) const
    {
        return _i == other._i and _j == other._j and _k == other._k;
    }
    
    bool operator!=(const GridIndex &other) const
    {
        return not operator==(other);
    }
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    bool isNull() const
    {
        return _k < 0;
    }
    
    static GridIndex null()
    {
        return GridIndex(0,0,-1);
    }
    
    qint32 i() const
    {
        return _i;
    }
    
    qint32 j() const
    {
        return _j;
    }
    
    qint32 k() const
    {
        return _k;
    }
    
private:
    qint32 _i, _j, _k;
};

/** This is a simple class that describes a 3D regular grid

    @author Christopher Woods
*/
class SIREVOL_EXPORT GridInfo
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const GridInfo&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, GridInfo&);

public:
    GridInfo();
    GridInfo(const AABox &dimensions, Length spacing);

    GridInfo(const GridInfo &other);

    ~GridInfo();
    
    GridInfo& operator=(const GridInfo &other);
    
    bool operator==(const GridInfo &other) const;
    bool operator!=(const GridInfo &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    GridIndex operator[](int i) const;

    int operator[](const GridIndex &idx) const;
    int operator[](const Vector &point) const;
    
    GridIndex at(int i) const;

    int at(const GridIndex &idx) const;
    int at(int i, int j, int k) const;
    int at(const Vector &point) const;
    
    GridIndex getitem(int i) const;
    
    AABox box(int i) const;
    AABox box(const GridIndex &idx) const;
    AABox box(const Vector &point) const;
    
    Vector point(int i) const;
    Vector point(const GridIndex &idx) const;
    Vector point(const Vector &point) const;
    
    qint32 dimX() const;
    qint32 dimY() const;
    qint32 dimZ() const;
    
    int nPoints() const;
    int count() const;
    int size() const;
    
    bool isEmpty() const;
    
    Length spacing() const;
    AABox dimensions() const;

    bool contains(const Vector &point) const;
    
    GridIndex closestIndexTo(const Vector &point) const;

    int gridToArrayIndex(int i, int j, int k) const;
    int gridToArrayIndex(const GridIndex &idx) const;
    int pointToArrayIndex(const Vector &point) const;

    void pointToGridCorners(const Vector &point, QVector<int> &indicies) const;
    void pointToGridCorners(const Vector &point, QVector<int> &indicies,
                            QVector<float> &weights) const;

    int pointToGridCorners(const MultiFloat &x, const MultiFloat &y, const MultiFloat &z,
                           QVector<MultiInt> &indicies) const;
    int pointToGridCorners(const MultiFloat &x, const MultiFloat &y, const MultiFloat &z,
                           QVector<MultiInt> &indicies, QVector<MultiFloat> &weights) const;

    GridIndex arrayToGridIndex(int i) const;
    GridIndex pointToGridIndex(const Vector &point) const;

    GridIndex indexOf(int i, const GridInfo &grid) const;
    GridIndex indexOf(const GridIndex &idx, const GridInfo &grid) const;

    QVector<float> redimension(const QVector<float> &values, const GridInfo &new_grid) const;

private:
    /** The origin of the grid */
    Vector grid_origin;
    
    /** The spacing of the grid */
    float grid_spacing;
    
    /** Inverse of the grid spacing */
    float inv_grid_spacing;
    
    /** The number of grid points along x, y, and z */
    qint32 dimx, dimy, dimz;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Get the grid index from the passed array index */
SIRE_ALWAYS_INLINE GridIndex GridInfo::arrayToGridIndex(int i) const
{
    if (i < 0 or i >= (dimx*dimy*dimz))
        return GridIndex::null();

    int ix = i / (dimy*dimz);
    i -= ix*dimy*dimz;

    int iy = i / dimz;
    i -= iy*dimz;

    return GridIndex(ix,iy,i);
}

/** Get the index into the 1D array for the point corresponding to grid 
    indicies i,j,k. */
SIRE_ALWAYS_INLINE int GridInfo::gridToArrayIndex(int i, int j, int k) const
{
    if (i < 0 or i >= dimx or j < 0 or j >= dimy or k < 0 or k >= dimz)
        return -1;
    else
        return i*(dimy*dimz) + j*dimz + k;
}

/** Get the index into the 1D array for the point corresponding to 
    the grid point with index 'idx'. Note that this returns '-1' if
    the passed GridIndex is null */
SIRE_ALWAYS_INLINE int GridInfo::gridToArrayIndex(const GridIndex &idx) const
{
    if (idx.isNull() or idx.i() >= dimx or idx.j() >= dimy or idx.k() >= dimz)
        return -1;
    else
        return idx.i()*(dimy*dimz) + idx.j()*dimz + idx.k();
}

/** Return the linear index of the point with grid index 'idx'. Note that this
    returns '-1' if the passed GridIndex is null */
SIRE_ALWAYS_INLINE int GridInfo::operator[](const GridIndex &idx) const
{
    return gridToArrayIndex(idx);
}

/** Return the grid index of the ith grid point */
SIRE_ALWAYS_INLINE GridIndex GridInfo::operator[](int i) const
{
    return arrayToGridIndex(i);
}

/** Return the linear array index of the box that contains point 'point'.
    Note that this will return -1 if the point lies outside the grid */
SIRE_ALWAYS_INLINE int GridInfo::operator[](const Vector &point) const
{
    return pointToArrayIndex(point);
}

/** Return the linear index of the point with grid index [i,j,k] */
SIRE_ALWAYS_INLINE int GridInfo::at(int i, int j, int k) const
{
    return gridToArrayIndex(i,j,k);
}

/** Return the linear index of the point with grid index 'idx'. Note that this
    returns '-1' if the passed GridIndex is null */
SIRE_ALWAYS_INLINE int GridInfo::at(const GridIndex &idx) const
{
    return operator[](idx);
}

/** Return the grid index of the ith grid point */
SIRE_ALWAYS_INLINE GridIndex GridInfo::at(int i) const
{
    return operator[](i);
}

/** Return the linear array index of the box that contains point 'point'.
    Note that this will return -1 if the point lies outside the grid */
SIRE_ALWAYS_INLINE int GridInfo::at(const Vector &point) const
{
    return pointToArrayIndex(point);
}

/** Return the number of grid points along the x dimension */
SIRE_ALWAYS_INLINE qint32 GridInfo::dimX() const
{
    return dimx;
}

/** Return the number of grid points along the y dimension */
SIRE_ALWAYS_INLINE qint32 GridInfo::dimY() const
{
    return dimy;
}

/** Return the number of grid points along the z dimension */
SIRE_ALWAYS_INLINE qint32 GridInfo::dimZ() const
{
    return dimz;
}

/** Return the total number of grid points */
SIRE_ALWAYS_INLINE int GridInfo::nPoints() const
{
    return dimx * dimy * dimz;
}

/** Return the total number of grid points */
SIRE_ALWAYS_INLINE int GridInfo::count() const
{
    return nPoints();
}

/** Return the total number of grid points */
SIRE_ALWAYS_INLINE int GridInfo::size() const
{
    return nPoints();
}

/** Return whether or not this grid is empty */
SIRE_ALWAYS_INLINE bool GridInfo::isEmpty() const
{
    return dimx == 0;
}

/** Return the grid spacing */
SIRE_ALWAYS_INLINE Length GridInfo::spacing() const
{
    return Length(grid_spacing);
}

/** Return the dimensions of the grid */
SIRE_ALWAYS_INLINE AABox GridInfo::dimensions() const
{
    return AABox::from( grid_origin, grid_origin + Vector((dimx-1)*grid_spacing,
                                                          (dimy-1)*grid_spacing,
                                                          (dimz-1)*grid_spacing) );
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireVol::GridIndex )
Q_DECLARE_METATYPE( SireVol::GridInfo )

Q_DECLARE_TYPEINFO( SireVol::GridIndex, Q_MOVABLE_TYPE );
Q_DECLARE_TYPEINFO( SireVol::GridInfo, Q_MOVABLE_TYPE );

SIRE_EXPOSE_CLASS( SireVol::GridIndex )
SIRE_EXPOSE_CLASS( SireVol::GridInfo )

SIRE_END_HEADER

#endif
