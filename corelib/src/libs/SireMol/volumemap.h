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

#ifndef SIREMOL_VOLUMEMAP_H
#define SIREMOL_VOLUMEMAP_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireVol/gridinfo.h"

#include "SireUnits/dimensions.h"

#include "moleculegroup.h"

#include <QList>
#include <QVector>

SIRE_BEGIN_HEADER

namespace SireMol
{
class VolumeMap;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::VolumeMap&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::VolumeMap&);

namespace SireMol
{

class Element;

using SireVol::GridInfo;
using SireUnits::Dimension::Length;
using SireBase::PropertyMap;

namespace detail{ class VolumeMapData; }

/** This class provides a volume map. This is a 3D regular grid,
    with the average occupancy at each grid point recorded. Grid
    points are considered occupied if they are covered by at least
    one atom
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT VolumeMap : public SireBase::ConcreteProperty<VolumeMap,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const VolumeMap&);
friend QDataStream& ::operator>>(QDataStream&, VolumeMap&);

public:
    enum MapType
    {
        AVERAGE = 1,
        MAXIMUM = 2,
        SUM = 3
    };

    enum FillType
    {
        POINT_ATOMS = 1,
        VDW_RADIUS = 2,
        COVALENT_RADIUS = 3,
        BOND_ORDER_RADIUS = 4
    };

    VolumeMap();
    VolumeMap(bool skip_light_atoms);
    
    VolumeMap(const Length &grid_spacing, bool skip_light_atoms=false);
    VolumeMap(MapType map_type, bool skip_light_atoms=false);
    VolumeMap(FillType fill_type, bool skip_light_atoms=false);
    
    VolumeMap(const Length &grid_spacing, MapType map_type, bool skip_light_atoms=false);
    VolumeMap(const Length &grid_spacing, FillType fill_type, bool skip_light_atoms=false);
    VolumeMap(FillType fill_type, MapType map_type, bool skip_light_atoms=false);
    VolumeMap(const Length &grid_spacing, FillType fill_type, MapType map_type,
              bool skip_light_atoms=false);
    
    VolumeMap(const VolumeMap &other);
    
    ~VolumeMap();
    
    VolumeMap& operator=(const VolumeMap &other);
    
    bool operator==(const VolumeMap &other) const;
    bool operator!=(const VolumeMap &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    bool isEmpty() const;
    
    int nGridPoints() const;
    int nMaxGridPoints() const;
    
    void setNMaxGridPoints(int nmax);
    
    Length gridSpacing() const;
    void setGridSpacing(Length grid_spacing);
    
    MapType mapType() const;
    void setMapType(MapType map_type);
    
    FillType fillType() const;
    void setFillType(FillType fill_type);
    
    void setSkipLightAtoms(bool on);
    bool skipLightAtoms() const;
    
    const GridInfo& gridInfo() const;
    
    qint64 nSamples() const;
    
    const QVector<float>& occupancy() const;
    
    void add(const MoleculeView &molecule, const PropertyMap &map = PropertyMap());
    void add(const MoleculeView &mol0, const MoleculeView &mol1,
             const PropertyMap &map = PropertyMap());
    void add(const QList<PartialMolecule> &molecules, const PropertyMap &map = PropertyMap());
    
    void add(const Molecules &molecules, const PropertyMap &map = PropertyMap());
    void add(const Molecules &mols0, const Molecules &mols1,
             const PropertyMap &map = PropertyMap());
    void add(const QList<Molecules> &molecules, const PropertyMap &map = PropertyMap());
    
    void add(const MoleculeGroup &molecules, const PropertyMap &map = PropertyMap());
    void add(const MoleculeGroup &mols0, const MoleculeGroup &mols1,
             const PropertyMap &map = PropertyMap());
    void add(const QList<MoleculeGroup> &molecules, const PropertyMap &map = PropertyMap());

    void add(const VolumeMap &other);
    void add(const GridInfo &gridinfo, const QVector<float> &values);
    
    void clear();

    bool isMasked() const;

    Length maskDistance() const;
    QVector<Vector> maskPoints() const;

    void setMaskWithinDistance(Length dist, const Vector &point,
                               bool clear_points = true);

    void setMaskWithinDistance(Length dist, const MoleculeView &molecule,
                               const PropertyMap &map = PropertyMap());

    void setMaskWithinDistance(Length dist, const MoleculeView &molecule,
                               bool clear_points, const PropertyMap &map = PropertyMap());

    void clearMask();

private:
    void redimensionGrid(Length new_spacing);
    void applyMask();

    void presize(const Molecules &molecules, const PropertyMap &map);
    AABox presize(const Vector &coords, const SireMol::Element &element, const AABox &box) const;

    void beginEvaluation();
    void evaluate(const MoleculeView &molecule, const PropertyMap &map);
    void evaluate(const Molecules &molecules, const PropertyMap &map);
    void evaluate(const GridInfo &info, const QVector<float> &vals, qint64 n);
    void cancelEvaluation();
    void endEvaluation();

    /** Information about the grid */
    GridInfo grid_info;
    
    /** The grid spacing desired by the user */
    Length grid_spacing;
    
    /** The method used to build the map */
    MapType map_type;
    
    /** The method used to assign atoms to points */
    FillType fill_type;
    
    /** The actual occupancy map */
    QVector<float> occ;
    
    /** Pointer to the data object used during grid updates */
    detail::VolumeMapData *d;
    
    /** The number of samples in this map */
    qint64 nsamples;
    
    /** The maximum number of grid points - this limits the memory
        used by this map */
    qint32 max_grid_points;

    /** The set of points used to provide a distance mask */
    QVector<Vector> mask_points;
    
    /** The distance used for the distance mask */
    float mask_dist;

    /** Whether or not to exclude light atoms from the map */
    bool skip_light_atoms;
};

}

Q_DECLARE_METATYPE( SireMol::VolumeMap )

SIRE_EXPOSE_CLASS( SireMol::VolumeMap )

SIRE_END_HEADER

#endif
