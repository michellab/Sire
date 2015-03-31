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

#include "gridinfo.h"

#include "SireMaths/vector.h"
#include "SireMaths/multifloat.h"
#include "SireMaths/multiint.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireVol;
using namespace SireMaths;
using namespace SireStream;

/////////
///////// Implementation of GridIndex
/////////

static const RegisterMetaType<GridIndex> r_index(NO_ROOT);

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const GridIndex &idx)
{
    writeHeader(ds, r_index, 1);
    ds << idx.i() << idx.j() << idx.k();
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, GridIndex &idx)
{
    VersionID v = readHeader(ds, r_index);
    
    if (v == 1)
    {
        qint32 i,j,k;
        ds >> i >> j >> k;
        idx = GridIndex(i,j,k);
    }
    else
        throw version_error( v, "1", r_index, CODELOC );
    
    return ds;
}

const char* GridIndex::typeName()
{
    return "SireVol::GridIndex";
}

const char* GridIndex::what() const
{
    return GridIndex::typeName();
}

QString GridIndex::toString() const
{
    if (isNull())
        return QObject::tr("GridIndex::null");
    else
        return QObject::tr("GridIndex( %1, %2, %3 )")
                        .arg(_i).arg(_j).arg(_k);
}

/////////
///////// Implementation of GridInfo
/////////

static const RegisterMetaType<GridInfo> r_info(NO_ROOT);

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const GridInfo &info)
{
    writeHeader(ds, r_info, 1);
    
    ds << info.grid_origin << info.grid_spacing
       << info.dimx << info.dimy << info.dimz;

    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, GridInfo &info)
{
    VersionID v = readHeader(ds, r_info);
    
    if (v == 1)
    {
        ds >> info.grid_origin >> info.grid_spacing
           >> info.dimx >> info.dimy >> info.dimz;
        
        info.inv_grid_spacing = 1.0f / info.grid_spacing;
    }
    else
        throw version_error(v, "1", r_info, CODELOC);
    
    return ds;
}

/** Constructor */
GridInfo::GridInfo()
         : grid_origin(-0.5), grid_spacing(1.0), inv_grid_spacing(1.0), dimx(0), dimy(0), dimz(0)
{}

/** Construct a grid of specified dimensions and spacing */
GridInfo::GridInfo(const AABox &dimensions, Length spacing)
         : grid_spacing(spacing.value()), dimx(0), dimy(0), dimz(0)
{
    if (grid_spacing <= 0)
    {
        throw SireError::invalid_arg( QObject::tr(
                "You cannot construct a grid with zero or negative grid spacing! "
                "(%1 A)").arg(grid_spacing), CODELOC );
    }
    
    inv_grid_spacing = 1.0f / grid_spacing;
    Vector nboxes = (2.0f * inv_grid_spacing) * dimensions.halfExtents();
    
    dimx = 1 + quint32( std::ceil(nboxes.x()) );
    dimy = 1 + quint32( std::ceil(nboxes.y()) );
    dimz = 1 + quint32( std::ceil(nboxes.z()) );
    
    grid_origin = dimensions.minCoords();
}

/** Copy constructor */
GridInfo::GridInfo(const GridInfo &other)
         : grid_origin(other.grid_origin), grid_spacing(other.grid_spacing),
           inv_grid_spacing(other.inv_grid_spacing),
           dimx(other.dimx), dimy(other.dimy), dimz(other.dimz)
{}

/** Destructor */
GridInfo::~GridInfo()
{}

/** Copy assignment operator */
GridInfo& GridInfo::operator=(const GridInfo &other)
{
    if (this != &other)
    {
        grid_origin = other.grid_origin;
        grid_spacing = other.grid_spacing;
        inv_grid_spacing = other.inv_grid_spacing;
        dimx = other.dimx;
        dimy = other.dimy;
        dimz = other.dimz;
    }
    
    return *this;
}

/** Comparison operator */
bool GridInfo::operator==(const GridInfo &other) const
{
    return grid_origin == other.grid_origin and grid_spacing == other.grid_spacing and
           dimx == other.dimx and dimy == other.dimy and dimz == other.dimz;
}

/** Comparison operator */
bool GridInfo::operator!=(const GridInfo &other) const
{
    return not operator==(other);
}

const char* GridInfo::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GridInfo>() );
}

const char* GridInfo::what() const
{
    return GridInfo::typeName();
}

QString GridInfo::toString() const
{
    return QObject::tr("GridInfo{ dimensions() = %1, spacing = %2 A, [%3, %4, %5] }")
                .arg(dimensions().toString())
                .arg(grid_spacing)
                .arg(dimx).arg(dimy).arg(dimz);
}

GridIndex GridInfo::getitem(int i) const
{
    return this->at(i);
}

/** Return the index of the grid box that contains the point 'point'. Note
    that this returns a null index if the point is not in the grid */
GridIndex GridInfo::pointToGridIndex(const Vector &point) const
{
    const qint32 i = qint32( (point.x() - grid_origin.x()) * inv_grid_spacing );
    const qint32 j = qint32( (point.y() - grid_origin.y()) * inv_grid_spacing );
    const qint32 k = qint32( (point.z() - grid_origin.z()) * inv_grid_spacing );
    
    if (i < 0 or i >= (dimx-1) or
        j < 0 or j >= (dimy-1) or
        k < 0 or k >= (dimz-1))
    {
        //the point is outside of the grid
        return GridIndex::null();
    }
    else
    {
        return GridIndex(i,j,k);
    }
}

/** Return array indicies of the eight grid points that are on the corners of the 
    box that contains the point 'point'. This returns eight '-1' values if the 
    point does not lie in the grid */
void GridInfo::pointToGridCorners(const Vector &point, QVector<int> &indicies) const
{
    indicies.resize(8);
    
    int *ia = indicies.data();
    
    const qint32 r = qint32( (point.x() - grid_origin.x()) * inv_grid_spacing );
    const qint32 s = qint32( (point.y() - grid_origin.y()) * inv_grid_spacing );
    const qint32 t = qint32( (point.z() - grid_origin.z()) * inv_grid_spacing );
    
    if (r < 0 or r >= (dimx-1) or
        s < 0 or s >= (dimy-1) or
        t < 0 or t >= (dimz-1))
    {
        for (int i=0; i<8; ++i)
        {
            ia[i] = -1;
        }
    }
    else
    {
        //i*(dimy*dimz) + j*dimz + k;
        
        const qint32 dimyz = dimz * dimy;
        
        ia[0] = (r  )*dimyz + (s  )*dimz + (t  );
        ia[1] = (r  )*dimyz + (s  )*dimz + (t+1);
        ia[2] = (r  )*dimyz + (s+1)*dimz + (t  );
        ia[3] = (r  )*dimyz + (s+1)*dimz + (t+1);
        ia[4] = (r+1)*dimyz + (s  )*dimz + (t  );
        ia[5] = (r+1)*dimyz + (s  )*dimz + (t+1);
        ia[6] = (r+1)*dimyz + (s+1)*dimz + (t  );
        ia[7] = (r+1)*dimyz + (s+1)*dimz + (t+1);
    }
}

/** Return array indicies of the eight grid points that are on the corners of the 
    box that contains the point 'point'. This returns eight '-1' values if the 
    point does not lie in the grid. This also returns the weights of the eight
    points, using tri-linear interpolation based on the distance between the 
    point and each corner of the box */
void GridInfo::pointToGridCorners(const Vector &point, QVector<int> &indicies,
                                  QVector<float> &weights) const
{
    indicies.resize(8);
    weights.resize(8);
    
    int *ia = indicies.data();
    float *wa = weights.data();
    
    float dx = (point.x() - grid_origin.x()) * inv_grid_spacing;
    float dy = (point.y() - grid_origin.y()) * inv_grid_spacing;
    float dz = (point.z() - grid_origin.z()) * inv_grid_spacing;
    
    const qint32 r = qint32(dx);
    const qint32 s = qint32(dy);
    const qint32 t = qint32(dz);
    
    dx -= r;
    dy -= s;
    dz -= t;
    
    const float one_minus_dx = 1.0f - dx;
    const float one_minus_dy = 1.0f - dy;
    const float one_minus_dz = 1.0f - dz;
    
    if (r < 0 or r >= (dimx-1) or
        s < 0 or s >= (dimy-1) or
        t < 0 or t >= (dimz-1))
    {
        for (int i=0; i<8; ++i)
        {
            ia[i] = -1;
            wa[i] = 0;
        }
    }
    else
    {
        //use tri-linear interpolation to get the weights
        //
        // This is described in 
        //
        // Davis, Madura and McCammon, Comp. Phys. Comm., 62, 187-197, 1991
        //
        // phi(x,y,z) = phi(i  ,j  ,k  )*(1-R)(1-S)(1-T) +
        //              phi(i+1,j  ,k  )*(  R)(1-S)(1-T) +
        //              phi(i  ,j+1,k  )*(1-R)(  S)(1-T) +
        //              phi(i  ,j  ,k+1)*(1-R)(1-S)(  T) +
        //              phi(i+1,j+1,k  )*(  R)(  S)(1-T) +
        //              phi(i+1,j  ,k+1)*(  R)(1-S)(  T) +
        //              phi(i  ,j+1,k+1)*(1-R)(  S)(  T) +
        //              phi(i+1,j+1,k+1)*(  R)(  S)(  T) +
        //
        // where R, S and T are the coordinates of the atom in 
        // fractional grid coordinates from the point (i,j,k), e.g.
        // (0,0,0) is (i,j,k) and (1,1,1) is (i+1,j+1,k+1)

        //i*(dimy*dimz) + j*dimz + k;
        
        const qint32 dimyz = dimz * dimy;
        
        ia[0] = (r  )*dimyz + (s  )*dimz + (t  );
        ia[1] = (r  )*dimyz + (s  )*dimz + (t+1);
        ia[2] = (r  )*dimyz + (s+1)*dimz + (t  );
        ia[3] = (r  )*dimyz + (s+1)*dimz + (t+1);
        ia[4] = (r+1)*dimyz + (s  )*dimz + (t  );
        ia[5] = (r+1)*dimyz + (s  )*dimz + (t+1);
        ia[6] = (r+1)*dimyz + (s+1)*dimz + (t  );
        ia[7] = (r+1)*dimyz + (s+1)*dimz + (t+1);
        
        wa[0] = (one_minus_dx) * (one_minus_dy) * (one_minus_dz);
        wa[1] = (one_minus_dx) * (one_minus_dy) * (    dz      );
        wa[2] = (one_minus_dx) * (    dy      ) * (one_minus_dz);
        wa[3] = (one_minus_dx) * (    dy      ) * (    dz      );
        wa[4] = (    dx      ) * (one_minus_dy) * (one_minus_dz);
        wa[5] = (    dx      ) * (one_minus_dy) * (    dz      );
        wa[6] = (    dx      ) * (    dy      ) * (one_minus_dz);
        wa[7] = (    dx      ) * (    dy      ) * (    dz      );
    }
}

/** Return array indicies of the eight grid points that are on the corners of the 
    box that contains the point 'point'. This returns eight '-1' values if the 
    point does not lie in the grid. The return value is the number of points
    that are in the box */
int GridInfo::pointToGridCorners(const MultiFloat &x, const MultiFloat &y,
                                 const MultiFloat &z, QVector<MultiInt> &indicies) const
{
    indicies.resize(8);
    
    MultiInt *ia = indicies.data();
    
    const MultiFloat ox( grid_origin.x() );
    const MultiFloat oy( grid_origin.y() );
    const MultiFloat oz( grid_origin.z() );
    const MultiFloat inv_spacing( inv_grid_spacing );
    
    const MultiInt r = (x - ox) * inv_spacing;
    const MultiInt s = (y - oy) * inv_spacing;
    const MultiInt t = (z - oz) * inv_spacing;
    
    const MultiInt one(1);
    
    //i*(dimy*dimz) + j*dimz + k;
    const MultiInt idz(dimz);
    const MultiInt idyz(dimz * dimy);
    
    ia[0] = (r    )*idyz + (s    )*idz + (t    );
    ia[1] = (r    )*idyz + (s    )*idz + (t+one);
    ia[2] = (r    )*idyz + (s+one)*idz + (t    );
    ia[3] = (r    )*idyz + (s+one)*idz + (t+one);
    ia[4] = (r+one)*idyz + (s    )*idz + (t    );
    ia[5] = (r+one)*idyz + (s    )*idz + (t+one);
    ia[6] = (r+one)*idyz + (s+one)*idz + (t    );
    ia[7] = (r+one)*idyz + (s+one)*idz + (t+one);

    //now check that the points are in the box
    int n_in_box = MultiFloat::count();
    
    for (int i=0; i<MultiInt::count(); ++i)
    {
        if (r[i] < 0 or r[i] >= (dimx-1) or
            s[i] < 0 or s[i] >= (dimy-1) or
            t[i] < 0 or t[i] >= (dimz-1))
        {
            for (int j=0; j<8; ++j)
            {
                ia[j].set(i, -1);
            }
            
            n_in_box -= 1;
        }
    }
    
    return n_in_box;
}

/** Return array indicies of the eight grid points that are on the corners of the 
    box that contains the point 'point'. This returns eight '-1' values if the 
    point does not lie in the grid. This also returns the weights of the eight
    points, using tri-linear interpolation based on the distance between the 
    point and each corner of the box */
int GridInfo::pointToGridCorners(const MultiFloat &x, const MultiFloat &y,
                                 const MultiFloat &z, QVector<MultiInt> &indicies,
                                 QVector<MultiFloat> &weights) const
{
    indicies.resize(8);
    weights.resize(8);
    
    MultiInt *ia = indicies.data();
    MultiFloat *wa = weights.data();
    
    const MultiFloat ox(grid_origin.x());
    const MultiFloat oy(grid_origin.y());
    const MultiFloat oz(grid_origin.z());
    
    const MultiFloat inv_spacing(inv_grid_spacing);
    
    MultiFloat dx = (x - ox) * inv_spacing;
    MultiFloat dy = (y - oy) * inv_spacing;
    MultiFloat dz = (z - oz) * inv_spacing;
    
    const MultiInt r = dx;
    const MultiInt s = dy;
    const MultiInt t = dz;
    
    dx -= r;
    dy -= s;
    dz -= t;

    const MultiFloat one_minus_dx = MULTIFLOAT_ONE - dx;
    const MultiFloat one_minus_dy = MULTIFLOAT_ONE - dy;
    const MultiFloat one_minus_dz = MULTIFLOAT_ONE - dz;

    //use tri-linear interpolation to get the weights
    //
    // This is described in 
    //
    // Davis, Madura and McCammon, Comp. Phys. Comm., 62, 187-197, 1991
    //
    // phi(x,y,z) = phi(i  ,j  ,k  )*(1-R)(1-S)(1-T) +
    //              phi(i+1,j  ,k  )*(  R)(1-S)(1-T) +
    //              phi(i  ,j+1,k  )*(1-R)(  S)(1-T) +
    //              phi(i  ,j  ,k+1)*(1-R)(1-S)(  T) +
    //              phi(i+1,j+1,k  )*(  R)(  S)(1-T) +
    //              phi(i+1,j  ,k+1)*(  R)(1-S)(  T) +
    //              phi(i  ,j+1,k+1)*(1-R)(  S)(  T) +
    //              phi(i+1,j+1,k+1)*(  R)(  S)(  T) +
    //
    // where R, S and T are the coordinates of the atom in 
    // fractional grid coordinates from the point (i,j,k), e.g.
    // (0,0,0) is (i,j,k) and (1,1,1) is (i+1,j+1,k+1)

    //i*(dimy*dimz) + j*dimz + k;
    
    const MultiInt one(1);
    
    //i*(dimy*dimz) + j*dimz + k;
    const MultiInt idz(dimz);
    const MultiInt idyz(dimz * dimy);
    
    ia[0] = (r    )*idyz + (s    )*idz + (t    );
    ia[1] = (r    )*idyz + (s    )*idz + (t+one);
    ia[2] = (r    )*idyz + (s+one)*idz + (t    );
    ia[3] = (r    )*idyz + (s+one)*idz + (t+one);
    ia[4] = (r+one)*idyz + (s    )*idz + (t    );
    ia[5] = (r+one)*idyz + (s    )*idz + (t+one);
    ia[6] = (r+one)*idyz + (s+one)*idz + (t    );
    ia[7] = (r+one)*idyz + (s+one)*idz + (t+one);
    
    wa[0] = (one_minus_dx) * (one_minus_dy) * (one_minus_dz);
    wa[1] = (one_minus_dx) * (one_minus_dy) * (    dz      );
    wa[2] = (one_minus_dx) * (    dy      ) * (one_minus_dz);
    wa[3] = (one_minus_dx) * (    dy      ) * (    dz      );
    wa[4] = (    dx      ) * (one_minus_dy) * (one_minus_dz);
    wa[5] = (    dx      ) * (one_minus_dy) * (    dz      );
    wa[6] = (    dx      ) * (    dy      ) * (one_minus_dz);
    wa[7] = (    dx      ) * (    dy      ) * (    dz      );
    
    //now check that the points are in the box
    int n_in_box = MultiFloat::count();
    
    for (int i=0; i<MultiInt::count(); ++i)
    {
        if (r[i] < 0 or r[i] >= (dimx-1) or
            s[i] < 0 or s[i] >= (dimy-1) or
            t[i] < 0 or t[i] >= (dimz-1))
        {
            for (int j=0; j<8; ++j)
            {
                ia[j].set(i, -1);
                wa[j].set(i, 0);
            }
            
            n_in_box -= 1;
        }
    }
    
    return n_in_box;
}

/** Return the array index of the grid box that contains the point 'point'. 
    Note that this returns -1 if the point is not in the grid */
int GridInfo::pointToArrayIndex(const Vector &point) const
{
    return this->gridToArrayIndex( pointToGridIndex(point) );
}

/** Return the AABox that encompasses the grid box at point 'idx'.
    This returns an empty box if the point is not in the grid */
AABox GridInfo::box(const GridIndex &idx) const
{
    if (idx.isNull() or idx.i() >= dimx or idx.j() >= dimy or idx.k() >= dimz)
    {
        return AABox();
    }
    else
    {
        return AABox( grid_origin + Vector( (0.5+idx.i()) * grid_spacing,
                                            (0.5+idx.j()) * grid_spacing,
                                            (0.5+idx.k()) * grid_spacing ),
                      Vector(0.5*grid_spacing) );
    }
}

/** Return the AABox that encompasses the grid box at point 'i'. This
    returns an empty box if there is no such point in the grid */
AABox GridInfo::box(int i) const
{
    return box( arrayToGridIndex(i) );
}

/** Return the AABox that encompasses the grid box that contains the
    point 'p'. Note that returns a null AABox if the point is not in the grid */
AABox GridInfo::box(const Vector &point) const
{
    return box( pointToGridIndex(point) );
}

/** Return the point at the bottom, left, back (lowest i, j, k indicies) of the
    box at grid index 'idx'. This returns a zero vector if the point is invalid. */
Vector GridInfo::point(const GridIndex &idx) const
{
    if (idx.isNull() or idx.i() >= dimx or idx.j() >= dimy or idx.k() >= dimz)
    {
        return Vector(0);
    }
    else
    {
        return Vector( grid_origin + Vector(idx.i() * grid_spacing,
                                            idx.j() * grid_spacing,
                                            idx.k() * grid_spacing ) );
    }
}

/** Return the point at the bottom, left, back (lowest i, j, k indicies) of the
    ith box. This returns a zero vector if the point is invalid. */
Vector GridInfo::point(int i) const
{
    return point( arrayToGridIndex(i) );
}

/** Return the point at the bottom, left, back (lowest i, j, k indicies) of the
    box containing the point 'point'. This returns a zero vector if the point is invalid. */
Vector GridInfo::point(const Vector &p) const
{
    return point( pointToGridIndex(p) );
}

/** Return whether or not this grid contains the point 'point' */
bool GridInfo::contains(const Vector &point) const
{
    if (this->isEmpty())
        return false;
    
    const Vector mincoords = dimensions().minCoords();
    const Vector maxcoords = dimensions().maxCoords();
    
    for (int i=0; i<3; ++i)
    {
        if (point[i] < mincoords[i] or point[i] > maxcoords[i])
            return false;
    }
    
    return true;
}

/** Return the grid index that is closest to the point 'point'. Note that
    if 'point' lies outside the grid, then the closest grid index will 
    still be returned (it may just lie a long way from the point!) */
GridIndex GridInfo::closestIndexTo(const Vector &point) const
{
    if (isEmpty())
        return GridIndex::null();
    
    float i = (point.x() - grid_origin.x()) * inv_grid_spacing;
    float j = (point.y() - grid_origin.y()) * inv_grid_spacing;
    float k = (point.z() - grid_origin.z()) * inv_grid_spacing;
    
    if (i < 0)
        i = 0;
    
    if (j < 0)
        j = 0;
    
    if (k < 0)
        k = 0;
    
    if (i >= dimx-1)
        i = dimx-1;
    
    if (j >= dimy-1)
        j = dimy-1;
    
    if (k >= dimz-1)
        k = dimz-1;

    //round the index to the nearest integer
    return GridIndex( qint32(i+0.5), qint32(j+0.5), qint32(k+0.5) );
}

/** Return the index in grid 'grid' of point 'idx' in this grid. This returns
    a null index if this index does not exist in either grid. Note that this
    returns the index of the closest grid point if the grids do not exactly 
    line up */
GridIndex GridInfo::indexOf(const GridIndex &idx, const GridInfo &grid) const
{
    if (idx.isNull())
        return GridIndex::null();
    
    const Vector point = this->point(idx);
    
    if (grid.contains(point))
        return grid.closestIndexTo(point);
    else
        return GridIndex::null();
}

/** Return the index in grid 'grid' of point 'i' in this grid. This returns
    a null index if this index does not exist in either grid. Note that this
    returns the index of the closest grid point if the grids do not exactly 
    line up */
GridIndex GridInfo::indexOf(int i, const GridInfo &grid) const
{
    return indexOf( this->arrayToGridIndex(i), grid );
}

/** Return the values 'values' that map to this grid, redimensioned to map to
    the grid 'new_grid' */
QVector<float> GridInfo::redimension(const QVector<float> &vals, const GridInfo &new_grid) const
{
    if (vals.count() != this->nPoints())
        throw SireError::incompatible_error( QObject::tr(
                "Cannot redimension the passed values as the number of grid points (%1) "
                "does not agree with the number of elements in the values array (%2).")
                    .arg(nPoints()).arg(vals.count()), CODELOC );
    
    if (this->operator==(new_grid))
        //no change in grid
        return vals;
    
    if (this->isEmpty() or new_grid.isEmpty())
        return QVector<float>();

    //create space to hold the expanded grid
    QVector<float> new_vals(new_grid.nPoints(), 0.0);
        
    for (int i=0; i<this->nPoints(); ++i)
    {
        int new_idx = new_grid.gridToArrayIndex( this->indexOf(i, new_grid) );
        
        if (new_idx >= 0 and new_idx < new_grid.nPoints())
        {
            new_vals.data()[new_idx] += vals.constData()[i];
        }
    }

    return new_vals;
}
