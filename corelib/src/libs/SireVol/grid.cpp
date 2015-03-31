/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "grid.h"

#include "SireMaths/rotate.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireVol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

//////////
////////// Implementation of Grid
//////////

static const RegisterMetaType<Grid> r_grid( MAGIC_ONLY, Grid::typeName() );

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const Grid &grid)
{
    writeHeader(ds, r_grid, 1);
    
    SharedDataStream sds(ds);
    
    sds << grid.grid_points << grid.grid_weights
        << static_cast<const Property&>(grid);
    
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, Grid &grid)
{
    VersionID v = readHeader(ds, r_grid);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> grid.grid_points >> grid.grid_weights
            >> static_cast<Property&>(grid);
            
        grid.aabox = AABox(grid.grid_points);
    }
    else
        throw version_error(v, "1", r_grid, CODELOC);
        
    return ds;
}

/** Base constructor */
Grid::Grid() : Property()
{}

/** Copy constructor */
Grid::Grid(const Grid &other) 
     : Property(other), grid_points(other.grid_points),
       grid_weights(other.grid_weights), aabox(other.aabox)
{}

/** Destructor */
Grid::~Grid()
{}

const char* Grid::typeName()
{
    return "SireVol::Grid";
}

/** Copy assignment operator */
Grid& Grid::operator=(const Grid &other)
{
    if (this != &other)
    {
        grid_points = other.grid_points;
        grid_weights = other.grid_weights;
        aabox = other.aabox;
        
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Grid::operator==(const Grid &other) const
{
    return grid_weights == other.grid_weights and
           grid_points == other.grid_points and
           Property::operator==(other);
}

/** Comparison operator */
bool Grid::operator!=(const Grid &other) const
{
    return not Grid::operator==(other);
}

/** Internal function used to set the grid */
void Grid::setGrid(const QVector<Vector> &gridpoints)
{
    grid_points = gridpoints;
    grid_weights = QVector<double>();
    aabox = AABox(grid_points);
}

/** Internal function used to set the grid points and grid weights */
void Grid::setGrid(const QVector<Vector> &gridpoints, const QVector<double> &weights)
{
    if (weights.isEmpty())
    {
        this->setGrid(gridpoints);
    }
    
    if (gridpoints.count() != weights.count())
        throw SireError::program_bug( QObject::tr(
                "It is a mistake to pass in a different number of grid point weights "
                "(%1) to grid point coordinates (%2)!")
                    .arg(gridpoints.count()).arg(weights.count()), CODELOC);
                    
    grid_points = gridpoints;
    grid_weights = weights;
    aabox = AABox(grid_points);
}

/** Return the array of grid points */
const QVector<Vector>& Grid::points() const
{
    return grid_points;
}

/** Return the array of grid weights - this will be empty
    if all of the points are equally weighted */
const QVector<double>& Grid::weights() const
{
    return grid_weights;
}

/** Return whether this grid is empty */
bool Grid::isEmpty() const
{
    return grid_points.isEmpty();
}

/** Return the minimum coordinates of the grid */
Vector Grid::minCoords() const
{
    return aabox.minCoords();
}

/** Return the maximum coordinates of the grid */
Vector Grid::maxCoords() const
{
    return aabox.maxCoords();
}

/** Return the center of the grid */
Vector Grid::center() const
{
    return aabox.center();
}

/** Return whether or not the grid points have weights
    (i.e. they are not equally weighted) */
bool Grid::hasWeights() const
{
    return not grid_weights.isEmpty();
}

/** Return the number of points in the grid */
int Grid::nPoints() const
{
    return grid_points.count();
}

/** Return the number of points in the grid */
int Grid::count() const
{
    return nPoints();
}

/** Return a raw pointer to the array of grid point coordinates */
const Vector* Grid::data() const
{
    return grid_points.constData();
}

/** Return a raw pointer to the array of grid point coordinates */
const Vector* Grid::constData() const
{
    return grid_points.constData();
}

/** Return a raw pointer to the array of grid point weights */
const double* Grid::weightData() const
{
    return grid_weights.constData();
}

/** Return a raw pointer to the array of grid point weights */
const double* Grid::constWeightData() const
{
    return grid_weights.constData();
}

Q_GLOBAL_STATIC( NullGrid, nullGrid )

const NullGrid& Grid::null()
{
    return *(nullGrid());
}

//////////
////////// Implementation of NullGrid
//////////

static const RegisterMetaType<NullGrid> r_nullgrid;

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const NullGrid &nullgrid)
{
    writeHeader(ds, r_nullgrid, 1);
    
    ds << static_cast<const Grid&>(nullgrid);
    
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, NullGrid &nullgrid)
{
    VersionID v = readHeader(ds, r_nullgrid);
    
    if (v == 1)
    {
        ds >> static_cast<Grid&>(nullgrid);
    }
    else
        throw version_error(v, "1", r_nullgrid, CODELOC);
        
    return ds;
}

/** Constructor */
NullGrid::NullGrid() : ConcreteProperty<NullGrid,Grid>()
{}

/** Copy constructor */
NullGrid::NullGrid(const NullGrid &other) : ConcreteProperty<NullGrid,Grid>(other)
{}

/** Destructor */
NullGrid::~NullGrid()
{}

/** Copy assigment operator */
NullGrid& NullGrid::operator=(const NullGrid &other)
{
    Grid::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullGrid::operator==(const NullGrid &other) const
{
    return Grid::operator==(other);
}

/** Comparison operator */
bool NullGrid::operator!=(const NullGrid &other) const
{
    return Grid::operator!=(other);
}

const char* NullGrid::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullGrid>() );
}

/** Return a copy of this grid that has been translated by 'delta' */
GridPtr NullGrid::translate(const Vector &delta) const
{
    return *this;
}

/** Return a copy of this grid that has been recentered to 'center' */
GridPtr NullGrid::recenter(const Vector &center) const
{
    return *this;
}

/** Return a copy of this grid that has been rotated using the passed rotation 
    quaternion about 'center' */
GridPtr NullGrid::rotate(const Quaternion &quat, const Vector &center) const
{
    return *this;
}

/** Return a copy of this grid that has been rotated using the passed rotation 
    matrix about 'center' */
GridPtr NullGrid::rotate(const Matrix &rotmat, const Vector &center) const
{
    return this->rotate( Quaternion(rotmat), center );
}

/** Return a copy of this grid that has been scaled uniformly by 'scalefactor' */
GridPtr NullGrid::scale(double scalefactor) const
{
    return *this;
}

//////////
////////// Implementation of RegularGrid
//////////

static const RegisterMetaType<RegularGrid> r_reggrid;

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const RegularGrid &grid)
{
    writeHeader(ds, r_reggrid, 1);
    
    ds << grid.basis_vectors << grid.grid_spacing
       << grid.dimx << grid.dimy << grid.dimz
       << static_cast<const Grid&>(grid);
       
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, RegularGrid &grid)
{
    VersionID v = readHeader(ds, r_reggrid);
    
    if (v == 1)
    {
        ds >> grid.basis_vectors >> grid.grid_spacing
           >> grid.dimx >> grid.dimy >> grid.dimz
           >> static_cast<Grid&>(grid);
    }
    else
        throw version_error( v, "1", r_reggrid, CODELOC );
        
    return ds;
}

/** Null constructor */
RegularGrid::RegularGrid()
            : ConcreteProperty<RegularGrid,Grid>(),
              basis_vectors( Matrix::identity() ),
              grid_spacing(0), dimx(0), dimy(0), dimz(0)
{}

/** Construct a regular grid that spans from 'min' to 'max' using the
    passed three orthoganol basis vectors, using a grid spacing of 'spacing' */
RegularGrid::RegularGrid(const Vector &min, const Vector &max,
                         const Quaternion &basis, Length spacing)
            : ConcreteProperty<RegularGrid,Grid>(),
              basis_vectors(basis.toMatrix()), grid_spacing(spacing)
{
    if (grid_spacing == Length(0))
    {
        this->operator=( RegularGrid() );
        return;
    }

    Vector mincoords = min;
    Vector maxcoords = max;
    mincoords.setMin(max);
    maxcoords.setMax(min);
    
    //work out the number of grid points in each dimension
    int nx = 1 + ((maxcoords.x() - mincoords.x()) / spacing.value());
    int ny = 1 + ((maxcoords.y() - mincoords.y()) / spacing.value());
    int nz = 1 + ((maxcoords.z() - mincoords.z()) / spacing.value());

    int n = nx * ny * nz;
    
    dimx = nx;
    dimy = ny;
    dimz = nz;
    
    if (nx > 1500 or ny > 1500 or nz > 1500 or n > (320*320*320))
        throw SireError::unavailable_resource( QObject::tr(
                    "You cannot create a grid of more than 320^3 points, "
                    "or with a single side having more than 1500 points, as this "
                    "would take up too much memory!"), CODELOC );
        
    Vector start = mincoords;
    
    QVector<Vector> grid_points(n);
    grid_points.squeeze();
    Vector *grid_point = grid_points.data();
    
    Vector delta_x = grid_spacing.value() * basis_vectors.column0();
    Vector delta_y = grid_spacing.value() * basis_vectors.column1();
    Vector delta_z = grid_spacing.value() * basis_vectors.column2();

    Vector current_x = start;
    
    for (int i=0; i<nx; ++i)
    {
        Vector current_y = current_x;
        
        for (int j=0; j<ny; ++j)
        {
            Vector current_z = current_y;
            
            for (int k=0; k<nz; ++k)
            {
                *grid_point = current_z;
                ++grid_point;
                current_z += delta_z;
            }
            
            current_y += delta_y;
        }
        
        current_x += delta_x;
    }
    
    this->setGrid(grid_points);
}

/** Construct a regular grid that spans from 'min' to 'max', with a grid
    spacing of 'spacing' */
RegularGrid::RegularGrid(const Vector &min, const Vector &max, Length spacing)
            : ConcreteProperty<RegularGrid,Grid>()
{
    this->operator=( RegularGrid(min, max, Quaternion(), spacing) );
}

/** Construct a regular grid that spans from 'min' to 'max' using the
    passed three orthoganol basis vectors, using a grid spacing of 'spacing' */
RegularGrid::RegularGrid(const Vector &min, const Vector &max, 
                         const Matrix &basis, Length spacing)
            : ConcreteProperty<RegularGrid,Grid>()
{
    this->operator=( RegularGrid(min, max, Quaternion(basis), spacing) );
}

/** Construct a regular grid that is centered at 'center', has a total of 
    'npoints' points which are spaced using a grid spacing of 'spacing',
    and arranged along the three orthoganol basis vectors supplied
    in 'basis' */
RegularGrid::RegularGrid(const Vector &center, const Quaternion &basis, int npoints,
                         Length spacing)
            : ConcreteProperty<RegularGrid,Grid>(),
              basis_vectors(basis.toMatrix()), grid_spacing(spacing)
{
    if (npoints > 320)
        throw SireError::unavailable_resource( QObject::tr(
                    "You cannot create a grid of more than 320^3 points, as this "
                    "would take up too much memory!"), CODELOC );

    if (npoints == 0 or grid_spacing == Length(0))
    {
        this->operator=( RegularGrid() );
        return;
    }
    else
    {
        double half_npoints = (npoints-1) / 2;
        
        Vector delta_x = grid_spacing * basis_vectors.column0();
        Vector delta_y = grid_spacing * basis_vectors.column1();
        Vector delta_z = grid_spacing * basis_vectors.column2();
        
        Vector start = center - half_npoints * delta_x 
                              - half_npoints * delta_y
                              - half_npoints * delta_z;
        
        int n = npoints * npoints * npoints;
        
        dimx = npoints;
        dimy = npoints;
        dimz = npoints;
        
        QVector<Vector> grid_points(n);
        grid_points.squeeze();
        Vector *grid_point = grid_points.data();
        
        Vector current_x = start;
        
        for (int i=0; i<npoints; ++i)
        {
            Vector current_y = current_x;
        
            for (int j=0; j<npoints; ++j)
            {
                Vector current_z = current_y;
            
                for (int k=0; k<npoints; ++k)
                {
                    *grid_point = current_z;
                    current_z += delta_z;
                    ++grid_point;
                }
                
                current_y += delta_y;
            }
            
            current_x += delta_x;
        }
        
        this->setGrid(grid_points);
    }
}

/** Construct a regular grid that is centered at 'center', has a total of 
    'npoints' points which are spaced using a grid spacing of 'spacing' */
RegularGrid::RegularGrid(const Vector &center, int npoints, Length spacing)
            : ConcreteProperty<RegularGrid,Grid>()
{
    this->operator=( RegularGrid(center, Quaternion(), npoints, spacing) );
}

/** Construct a regular grid that is centered at 'center', has a total of 
    'npoints' points which are spaced using a grid spacing of 'spacing',
    and arranged along the three orthoganol basis vectors supplied
    in 'basis' */
RegularGrid::RegularGrid(const Vector &center, const Matrix &basis, int npoints,
                         Length spacing)
            : ConcreteProperty<RegularGrid,Grid>()
{
    this->operator=( RegularGrid(center, Quaternion(basis), npoints, spacing) );
}

/** Copy constructor */
RegularGrid::RegularGrid(const RegularGrid &other)
            : ConcreteProperty<RegularGrid,Grid>(other),
              basis_vectors(other.basis_vectors), grid_spacing(other.grid_spacing),
              dimx(other.dimx), dimy(other.dimy), dimz(other.dimz)
{}

/** Destructor */
RegularGrid::~RegularGrid()
{}

/** Copy assignment operator */
RegularGrid& RegularGrid::operator=(const RegularGrid &other)
{
    if (this != &other)
    {
        basis_vectors = other.basis_vectors;
        grid_spacing = other.grid_spacing;
        dimx = other.dimx;
        dimy = other.dimy;
        dimz = other.dimz;
        Grid::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool RegularGrid::operator==(const RegularGrid &other) const
{
    return basis_vectors == other.basis_vectors and
           grid_spacing == other.grid_spacing and
           dimx == other.dimx and dimy == other.dimy and
           dimz == other.dimz and Grid::operator==(other);
}

/** Comparison operator */
bool RegularGrid::operator!=(const RegularGrid &other) const
{
    return not RegularGrid::operator==(other);
}

/** Return a string representation of the grid */
QString RegularGrid::toString() const
{
    return QObject::tr("RegularGrid( from %1 to %2, spacing %3 A )")
                .arg(minCoords().toString(), maxCoords().toString())
                .arg(grid_spacing.to(angstrom));
}

/** Return the basis vectors for the grid */
const Matrix& RegularGrid::basis() const
{
    return basis_vectors;
}

/** Return the grid spacing (this is a cubic grid) */
Length RegularGrid::gridSpacing() const
{
    return grid_spacing;
}

/** Return the number of points in the x dimension */
int RegularGrid::dimX() const
{
    return dimx;
}

/** Return the number of points in the y dimension */
int RegularGrid::dimY() const
{
    return dimy;
}

/** Return the number of points in the z dimension */
int RegularGrid::dimZ() const
{
    return dimz;
}

/** Return a copy of this grid that has been translated by 'delta' */
GridPtr RegularGrid::translate(const Vector &delta) const
{
    if (this->isEmpty() or delta.isZero())
        return *this;

    RegularGrid grid(*this);

    QVector<Vector> grid_points = this->points();
    
    for (QVector<Vector>::iterator it = grid_points.begin();
         it != grid_points.end();
         ++it)
    {
        *it += delta;
    }

    grid.setGrid(grid_points);

    return grid;
}

/** Return a copy of this grid that has been recentered to 'center' */
GridPtr RegularGrid::recenter(const Vector &cent) const
{
    return this->translate( cent - this->center() );
}

/** Return a copy of this grid that has been rotated using the passed rotation 
    quaternion about 'center' */
GridPtr RegularGrid::rotate(const Quaternion &quat, const Vector &center) const
{
    if (this->isEmpty())
        return *this;

    RegularGrid grid(*this);

    QVector<Vector> grid_points = this->points();
    
    Matrix rotmat = quat.toMatrix();
    
    for (QVector<Vector>::iterator it = grid_points.begin();
         it != grid_points.end();
         ++it)
    {
        *it = SireMaths::rotate(*it, rotmat, center);
    }

    grid.setGrid(grid_points);
    grid.basis_vectors *= rotmat;

    return grid;
}

/** Return a copy of this grid that has been rotated using the passed rotation 
    matrix about 'center' */
GridPtr RegularGrid::rotate(const Matrix &rotmat, const Vector &center) const
{
    return this->rotate( Quaternion(rotmat), center );
}

/** Return a copy of this grid that has been scaled uniformly by 'scalefactor' */
GridPtr RegularGrid::scale(double scalefactor) const
{
    if (this->isEmpty() or scalefactor == 1)
        return *this;

    RegularGrid grid(*this);

    QVector<Vector> grid_points = this->points();
    
    Vector center = this->center();
    
    for (QVector<Vector>::iterator it = grid_points.begin();
         it != grid_points.end();
         ++it)
    {
        *it = center + scalefactor * (*it - center);
    }

    grid.setGrid(grid_points);
    grid.grid_spacing = Length(grid_spacing.value() * scalefactor);

    return grid;
}
