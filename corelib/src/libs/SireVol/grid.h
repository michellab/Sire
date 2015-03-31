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

#ifndef SIREVOL_GRID_H
#define SIREVOL_GRID_H

#include "SireBase/property.h"

#include "SireMaths/vector.h"
#include "SireMaths/matrix.h"
#include "SireMaths/quaternion.h"

#include "SireUnits/dimensions.h"

#include "aabox.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class Grid;
class NullGrid;

class RegularGrid;
}

QDataStream& operator<<(QDataStream&, const SireVol::Grid&);
QDataStream& operator>>(QDataStream&, SireVol::Grid&);

QDataStream& operator<<(QDataStream&, const SireVol::NullGrid&);
QDataStream& operator>>(QDataStream&, SireVol::NullGrid&);

QDataStream& operator<<(QDataStream&, const SireVol::RegularGrid&);
QDataStream& operator>>(QDataStream&, SireVol::RegularGrid&);

namespace SireVol
{

typedef SireBase::PropPtr<Grid> GridPtr;

using SireMaths::Vector;
using SireMaths::Quaternion;
using SireMaths::Matrix;

/** This is the base class of all grids. A grid is a set of points
    in space which are laid out in some manner (e.g. regularly spaced,
    or randomly distributed). 
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT Grid : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const Grid&);
friend QDataStream& ::operator>>(QDataStream&, Grid&);

public:
    Grid();
    Grid(const Grid &other);
    
    ~Grid();
    
    static const char* typeName();
    
    virtual Grid* clone() const=0;
    
    const QVector<Vector>& points() const;
    const QVector<double>& weights() const;

    bool isEmpty() const;

    Vector minCoords() const;
    Vector maxCoords() const;
    Vector center() const;
    
    const AABox& aaBox() const;
    
    bool hasWeights() const;
    
    int nPoints() const;
    
    int count() const;
    
    const Vector* data() const;
    const Vector* constData() const;
    
    const double* weightData() const;
    const double* constWeightData() const;
    
    static const NullGrid& null();
    
    virtual GridPtr translate(const Vector &delta) const=0;

    virtual GridPtr recenter(const Vector &center) const=0;
    
    virtual GridPtr rotate(const Matrix &rotmat, 
                           const Vector &center = Vector(0)) const=0;
    virtual GridPtr rotate(const Quaternion &quat,
                           const Vector &center = Vector(0)) const=0;
    
    virtual GridPtr scale(double scalefactor) const=0;
    
protected:
    Grid& operator=(const Grid &other);
    
    bool operator==(const Grid &other) const;
    bool operator!=(const Grid &other) const;

    void setGrid(const QVector<Vector> &gridpoints);
    void setGrid(const QVector<Vector> &gridpoints,
                 const QVector<double> &weights);

private:
    /** The locations of all of the grid points */
    QVector<Vector> grid_points;
    
    /** The weights of all of the points - this is empty
        if all of the points are equally weighted */
    QVector<double> grid_weights;
    
    /** The axis-aligned box that encompasses all of the grid points */
    AABox aabox;
};

/** This is the null (empty) grid */
class SIREVOL_EXPORT NullGrid : public SireBase::ConcreteProperty<NullGrid,Grid>
{

friend QDataStream& ::operator<<(QDataStream&, const NullGrid&);
friend QDataStream& ::operator>>(QDataStream&, NullGrid&);

public:
    NullGrid();
    NullGrid(const NullGrid &other);
    
    ~NullGrid();
    
    NullGrid& operator=(const NullGrid &other);
    
    bool operator==(const NullGrid &other) const;
    bool operator!=(const NullGrid &other) const;
    
    static const char* typeName();
    
    GridPtr translate(const Vector &delta) const;
    
    GridPtr recenter(const Vector &center) const;
    
    GridPtr rotate(const Matrix &rotmat, const Vector &center = Vector(0)) const;
    GridPtr rotate(const Quaternion &quat, const Vector &center = Vector(0)) const;
    
    GridPtr scale(double scalefactor) const;
};

/** This is a regular grid, containing points which are equally spaced
    along the three orthoganol coordinates
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT RegularGrid : public SireBase::ConcreteProperty<RegularGrid,Grid>
{

friend QDataStream& ::operator<<(QDataStream&, const RegularGrid&);
friend QDataStream& ::operator>>(QDataStream&, RegularGrid&);

public:
    RegularGrid();
    RegularGrid(const Vector &min, const Vector &max, 
                SireUnits::Dimension::Length gridsize);
    RegularGrid(const Vector &min, const Vector &max, 
                const Matrix &basis, SireUnits::Dimension::Length gridsize);
    RegularGrid(const Vector &min, const Vector &max,
                const Quaternion &basis, SireUnits::Dimension::Length gridsize);

    RegularGrid(const Vector &center, int npoints,
                SireUnits::Dimension::Length gridsize);
    RegularGrid(const Vector &center, const Matrix &basis, int npoints,
                SireUnits::Dimension::Length gridsize);
    RegularGrid(const Vector &center, const Quaternion &basis, int npoints,
                SireUnits::Dimension::Length gridsize);

    RegularGrid(const RegularGrid &other);
    
    ~RegularGrid();
    
    RegularGrid& operator=(const RegularGrid &other);
    
    bool operator==(const RegularGrid &other) const;
    bool operator!=(const RegularGrid &other) const;

    QString toString() const;
    
    const Matrix& basis() const;
    
    int dimX() const;
    int dimY() const;
    int dimZ() const;
    
    SireUnits::Dimension::Length gridSpacing() const;
    
    GridPtr translate(const Vector &delta) const;
    
    GridPtr recenter(const Vector &center) const;
    
    GridPtr rotate(const Matrix &rotmat, const Vector &center = Vector(0)) const;
    GridPtr rotate(const Quaternion &quat, const Vector &center = Vector(0)) const;
    
    GridPtr scale(double scalefactor) const;
    
private:
    /** The basis vectors for this grid */
    Matrix basis_vectors;
    
    /** The distance between nearest neighbour grid points */
    SireUnits::Dimension::Length grid_spacing;
    
    /** The number of points in the x, y and z dimensions */
    quint32 dimx, dimy, dimz;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the axis-aligned box that just encompasses the grid points */
inline const AABox& Grid::aaBox() const
{
    return aabox;
}

#endif

}

Q_DECLARE_METATYPE( SireVol::NullGrid )
Q_DECLARE_METATYPE( SireVol::RegularGrid )

SIRE_EXPOSE_CLASS( SireVol::Grid )
SIRE_EXPOSE_CLASS( SireVol::NullGrid )
SIRE_EXPOSE_CLASS( SireVol::RegularGrid )

SIRE_EXPOSE_PROPERTY( SireVol::GridPtr, SireVol::Grid )

SIRE_END_HEADER

#endif
