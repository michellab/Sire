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

#include "volumemap.h"

#include "molecule.h"
#include "partialmolecule.h"
#include "atomcoords.h"
#include "atomelements.h"

#include "SireVol/cartesian.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>

using namespace SireMol;
using namespace SireBase;
using namespace SireVol;
using namespace SireUnits;
using namespace SireStream;

static const RegisterMetaType<VolumeMap> r_volmap;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const VolumeMap &volmap)
{
    writeHeader(ds, r_volmap, 2);
    
    SharedDataStream sds(ds);
    
    sds << volmap.grid_info << volmap.grid_spacing.to(angstrom)
        << qint32(volmap.map_type) << qint32(volmap.fill_type)
        << volmap.occ << volmap.nsamples << volmap.max_grid_points
        << volmap.mask_points << volmap.mask_dist
        << volmap.skip_light_atoms
        << static_cast<const Property&>(volmap);
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, VolumeMap &volmap)
{
    VersionID v = readHeader(ds, r_volmap);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        double grid_spacing;
        qint32 map_type, fill_type;
        
        sds >> volmap.grid_info >> grid_spacing
            >> map_type >> fill_type
            >> volmap.occ >> volmap.nsamples >> volmap.max_grid_points
            >> volmap.mask_points >> volmap.mask_dist
            >> volmap.skip_light_atoms
            >> static_cast<Property&>(volmap);
        
        volmap.grid_spacing = grid_spacing * angstrom;
        volmap.map_type = VolumeMap::MapType(map_type);
        volmap.fill_type = VolumeMap::FillType(fill_type);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        double grid_spacing;
        qint32 map_type, fill_type;
        
        volmap.mask_dist = 0;
        volmap.mask_points.clear();
        
        sds >> volmap.grid_info >> grid_spacing
            >> map_type >> fill_type
            >> volmap.occ >> volmap.nsamples >> volmap.max_grid_points
            >> volmap.skip_light_atoms
            >> static_cast<Property&>(volmap);
        
        volmap.grid_spacing = grid_spacing * angstrom;
        volmap.map_type = VolumeMap::MapType(map_type);
        volmap.fill_type = VolumeMap::FillType(fill_type);
    }
    else
        throw version_error(v, "1,2", r_volmap, CODELOC);
    
    return ds;
}

namespace SireMol
{
    namespace detail
    {
        class VolumeMapData
        {
        public:
            VolumeMapData()
            {}
            
            ~VolumeMapData()
            {}
            
            void addToGrid(const Vector &coords, const Element &element);
            
            GridInfo grid_info;
            Length grid_spacing;
            QVector<float> occ;
            
            VolumeMap::MapType map_type;
            VolumeMap::FillType fill_type;
            
            QVector<Vector> mask_points;
            float mask_dist;
            
            qint32 max_points;
            
            QElapsedTimer t;
        };
    }
}

static const Length default_grid_spacing = 0.5 * angstrom;
static const VolumeMap::MapType default_map_type = VolumeMap::AVERAGE;
static const VolumeMap::FillType default_fill_type = VolumeMap::VDW_RADIUS;
static const qint32 default_max_grid_points = 256*256*256;  // 64 MB!!!, but 128^3 A^3 at default

static const Length min_grid_spacing = 0.01 * angstrom;

/** Default constructor */
VolumeMap::VolumeMap() : ConcreteProperty<VolumeMap,Property>(),
                         grid_spacing(default_grid_spacing),
                         map_type(default_map_type),
                         fill_type(default_fill_type),
                         d(0), nsamples(0), max_grid_points(default_max_grid_points),
                         mask_dist(0),
                         skip_light_atoms(false)
{}

/** Default constructor */
VolumeMap::VolumeMap(bool skip) : ConcreteProperty<VolumeMap,Property>(),
                                  grid_spacing(default_grid_spacing),
                                  map_type(default_map_type),
                                  fill_type(default_fill_type),
                                  d(0), nsamples(0), max_grid_points(default_max_grid_points),
                                  mask_dist(0),
                                  skip_light_atoms(skip)
{}

/** Construct, specifying the grid spacing */
VolumeMap::VolumeMap(const Length &spacing, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(spacing),
            map_type(default_map_type),
            fill_type(default_fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{
    if (grid_spacing.value() < min_grid_spacing)
    {
        grid_spacing = min_grid_spacing;
    }
}

/** Construct, specifying the type of map to build */
VolumeMap::VolumeMap(MapType map_type, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(default_grid_spacing),
            map_type(map_type),
            fill_type(default_fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{}

/** Construct, specifying the method used to assign atoms to grid points */
VolumeMap::VolumeMap(FillType fill_type, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(default_grid_spacing),
            map_type(default_map_type),
            fill_type(fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{}

/** Construct, specifying the grid spacing and type of map to build */
VolumeMap::VolumeMap(const Length &spacing, MapType map_type, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(spacing),
            map_type(map_type),
            fill_type(default_fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{
    if (grid_spacing.value() < min_grid_spacing)
    {
        grid_spacing = min_grid_spacing;
    }
}

/** Construct, specifying the grid spacing and method of assigning atoms to points */
VolumeMap::VolumeMap(const Length &spacing, FillType fill_type, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(spacing),
            map_type(default_map_type),
            fill_type(fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{
    if (grid_spacing.value() < min_grid_spacing)
    {
        grid_spacing = min_grid_spacing;
    }
}

/** Construct, specifying the type of map and method of assigning atoms to points */
VolumeMap::VolumeMap(FillType fill_type, MapType map_type, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(default_grid_spacing),
            map_type(map_type),
            fill_type(fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{}

/** Construct, specifying the grid spacing, map type and method of assigning atoms to points */
VolumeMap::VolumeMap(const Length &spacing, FillType fill_type, MapType map_type, bool skip)
          : ConcreteProperty<VolumeMap,Property>(),
            grid_spacing(spacing),
            map_type(map_type),
            fill_type(fill_type),
            d(0), nsamples(0), max_grid_points(default_max_grid_points),
            mask_dist(0),
            skip_light_atoms(skip)
{
    if (grid_spacing.value() < min_grid_spacing)
    {
        grid_spacing = min_grid_spacing;
    }
}

/** Copy constructor */
VolumeMap::VolumeMap(const VolumeMap &other)
          : ConcreteProperty<VolumeMap,Property>(other),
            grid_info(other.grid_info), grid_spacing(other.grid_spacing),
            map_type(other.map_type), fill_type(other.fill_type),
            occ(other.occ), d(0), nsamples(other.nsamples), max_grid_points(other.max_grid_points),
            mask_points(other.mask_points), mask_dist(other.mask_dist),
            skip_light_atoms(other.skip_light_atoms)
{}

/** Destructor */
VolumeMap::~VolumeMap()
{
    delete d;
}

/** Copy assignment operator */
VolumeMap& VolumeMap::operator=(const VolumeMap &other)
{
    if (this != &other)
    {
        grid_info = other.grid_info;
        grid_spacing = other.grid_spacing;
        map_type = other.map_type;
        fill_type = other.fill_type;
        occ = other.occ;
        nsamples = other.nsamples;
        max_grid_points = other.max_grid_points;
        mask_points = other.mask_points;
        mask_dist = other.mask_dist;
        skip_light_atoms = other.skip_light_atoms;
        
        delete d;
        d = 0;
    }
    
    return *this;
}

/** Comparison operator */
bool VolumeMap::operator==(const VolumeMap &other) const
{
    return this == &other or
           (grid_info == other.grid_info and grid_spacing == other.grid_spacing and
            map_type == other.map_type and fill_type == other.fill_type and
            occ == other.occ and nsamples == other.nsamples and
            max_grid_points == other.max_grid_points and
            mask_points == other.mask_points and
            mask_dist == other.mask_dist and
            skip_light_atoms == other.skip_light_atoms);
}

/** Comparison operator */
bool VolumeMap::operator!=(const VolumeMap &other) const
{
    return not operator==(other);
}

const char* VolumeMap::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VolumeMap>() );
}

const char* VolumeMap::what() const
{
    return VolumeMap::typeName();
}

/** Return whether or not this map is empty */
bool VolumeMap::isEmpty() const
{
    return nsamples == 0;
}

QString VolumeMap::toString() const
{
    if (isEmpty())
        return QObject::tr("VolumeMap( gridSpacing() == %1 A )").arg(grid_spacing.to(angstrom));
    else
        return QObject::tr("VolumeMap( gridInfo() == %1, nSamples() == %2, "
                           "skipLightAtoms() == %3 )")
                    .arg(gridInfo().toString()).arg(nsamples).arg(skipLightAtoms());
}

/** Return the number of grid points in the grid */
int VolumeMap::nGridPoints() const
{
    return grid_info.nPoints();
}

/** Return the maximum number of grid points available to this map */
int VolumeMap::nMaxGridPoints() const
{
    return max_grid_points;
}

/** Set the maximum number of grid points available to this map. Note that
    if this is less than the current number of points, then it prevents
    the map from growing (but doesn't shrink the map) */
void VolumeMap::setNMaxGridPoints(int nmax)
{
    max_grid_points = nmax;
}

/** Return the grid spacing */
Length VolumeMap::gridSpacing() const
{
    return grid_spacing;
}

/** Set whether or not to exclude light atoms from the map */
void VolumeMap::setSkipLightAtoms(bool skip)
{
    skip_light_atoms = skip;
}

/** Return whether or not light atoms are excluded from the map */
bool VolumeMap::skipLightAtoms() const
{
    return skip_light_atoms;
}

/** Clear the current grid */
void VolumeMap::clear()
{
    nsamples = 0;
    occ = QVector<float>();
    grid_info = GridInfo();
    delete d;
    d = 0;
}

/** Internal function called to change the grid spacing of the current grid */
void VolumeMap::redimensionGrid(Length new_spacing)
{
    GridInfo new_grid_info( grid_info.dimensions(), new_spacing );
    
    if (new_grid_info.nPoints() > max_grid_points)
        throw SireError::unavailable_resource( QObject::tr(
                "Cannot redimension this grid to a spacing %1 A, as doing so would "
                "increase the number of required grid points to %2, which is more than "
                "the maximum allowed. Either choose a larger grid spacing or increase "
                "the maximum number of allowable grid points.")
                    .arg(new_spacing.value())
                    .arg(new_grid_info.nPoints())
                    .arg(max_grid_points), CODELOC );
    
    if (not occ.isEmpty())
    {
        occ = grid_info.redimension(occ, new_grid_info);
    }
    
    grid_info = new_grid_info;
    grid_spacing = new_spacing;
}

/** Set the desired grid spacing. If this is not the same as the current
    grid, then the current map is re-mapped onto the new grid */
void VolumeMap::setGridSpacing(Length new_spacing)
{
    if (new_spacing < min_grid_spacing)
        new_spacing = min_grid_spacing;
    
    if (grid_spacing != new_spacing)
    {
        redimensionGrid(new_spacing);
    }
    
    grid_spacing = new_spacing;
}

/** Return the type of map */
VolumeMap::MapType VolumeMap::mapType() const
{
    return map_type;
}

/** Set the type of map. Supported types are;

    AVERAGE - accumulates the average occupancy of each point
    MAXIMUM - accumulates the maximum occupancy of each point
    SUM     - accumulates the summed occupancy of each point

    Default is AVERAGE

    Note that changing the map type will clear the current grid
*/
void VolumeMap::setMapType(MapType new_type)
{
    if (new_type != map_type)
        clear();
    
    map_type = new_type;
}

/** Return the method used to assign atoms to grid points */
VolumeMap::FillType VolumeMap::fillType() const
{
    return fill_type;
}

/** Set the method to assign atoms to grid points. Supported methods are;

    POINT_ATOMS - atoms are only assigned to the grid point closest to their center
    VDW_RADIUS  - atoms are assigned to all points underneath their VDW radius
    COVALENT_RADIUS - atoms are assigned to all points underneath their covalent radius
    BOND_ORDER_RADIUS - atoms are assigned to all points underneath their bond order radius
    
    Default is VDW_RADIUS
    
    Note that changing the fill type will clear the current grid
*/
void VolumeMap::setFillType(FillType new_type)
{
    if (new_type != fill_type)
        clear();
    
    fill_type = new_type;
}

/** Return information about the grid. Note that the grid will grow automatically
    to cover atoms as they are added to the map */
const GridInfo& VolumeMap::gridInfo() const
{
    return grid_info;
}

/** Return the current occupancy map. This array of values should be read in conjunction
    with the current GridInfo */
const QVector<float>& VolumeMap::occupancy() const
{
    return occ;
}

/** Return the number of samples used to create this map */
qint64 VolumeMap::nSamples() const
{
    return nsamples;
}

AABox VolumeMap::presize(const Vector &coords, const Element &element, const AABox &box) const
{
    float rad;
    
    switch (fill_type)
    {
        case VolumeMap::VDW_RADIUS:
            rad = element.vdwRadius();
            break;
        case VolumeMap::COVALENT_RADIUS:
            rad = element.covalentRadius();
            break;
        case VolumeMap::BOND_ORDER_RADIUS:
            rad = element.bondOrderRadius();
            break;
        case VolumeMap::POINT_ATOMS:
        default:
            rad = 0;
    }
    
    //ensure that the grid contains the atom
    const AABox atombox = AABox(coords, Vector(rad+2.5));

    if (box.isNull())
        return atombox;
    
    else if (box.contains(atombox))
        return box;
    
    else
    {
        return box + atombox;
    }
}

/** Internal function called to presize the grid (to prevent lots of memory
    reallocation when adding in the first block of molecules) */
void VolumeMap::presize(const Molecules &molecules, const PropertyMap &map)
{
    if (not occ.isEmpty())
        return;

    AABox box;
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        const MoleculeView &molecule = it.value();
    
        //extract the coordinates and element properties from the molecule
        const AtomCoords &coords = molecule.data().property( map["coordinates"] ).asA<AtomCoords>();
        const AtomElements &elems = molecule.data().property( map["element"] ).asA<AtomElements>();
        
        if (molecule.selectedAll())
        {
            //loop over all atoms and add them to the grid
            for (int i=0; i<coords.nCutGroups(); ++i)
            {
                const Vector *ca = coords.data(CGIdx(i));
                const Element *ea = elems.data(CGIdx(i));
                
                for (int j=0; j<coords.nAtoms(CGIdx(i)); ++j)
                {
                    if ((not skipLightAtoms()) or ea[j].nProtons() >= 6)
                    {
                        box = presize(ca[j], ea[j], box);
                    }
                }
            }
        }
        else
        {
            //not selected all atoms...
            const AtomSelection selected_atoms = molecule.selection();
            
            for (int i=0; i<coords.nCutGroups(); ++i)
            {
                const Vector *ca = coords.data(CGIdx(i));
                const Element *ea = elems.data(CGIdx(i));

                if (selected_atoms.selectedAll(CGIdx(i)))
                {
                    for (int j=0; j<coords.nAtoms(CGIdx(i)); ++j)
                    {
                        if ((not skipLightAtoms()) or ea[j].nProtons() >= 6)
                        {
                            box = presize(ca[j], ea[j], box);
                        }
                    }
                }
                else
                {
                    foreach (Index j, selected_atoms.selectedAtoms(CGIdx(i)))
                    {
                        if ((not skipLightAtoms()) or ea[j.value()].nProtons() >= 6)
                        {
                            box = presize(ca[j], ea[j], box);
                        }
                    }
                }
            }
        }
    }
}

/** Begin a new evaluation */
void VolumeMap::beginEvaluation()
{
    delete d;
    d = 0;
    
    d = new detail::VolumeMapData();
    
    d->grid_info = grid_info;
    d->occ = QVector<float>(grid_info.nPoints(), 0.0);
    d->grid_spacing = grid_spacing;
    d->map_type = map_type;
    d->fill_type = fill_type;
    d->max_points = max_grid_points;
    d->mask_points = mask_points;
    d->mask_dist = mask_dist;
    //d->t.start();
}

void SireMol::detail::VolumeMapData::addToGrid(const Vector &coords, const Element &element)
{
    float rad;
    
    switch (fill_type)
    {
        case VolumeMap::VDW_RADIUS:
            rad = element.vdwRadius();
            break;
        case VolumeMap::COVALENT_RADIUS:
            rad = element.covalentRadius();
            break;
        case VolumeMap::BOND_ORDER_RADIUS:
            rad = element.bondOrderRadius();
            break;
        case VolumeMap::POINT_ATOMS:
        default:
            rad = 0;
    }

    const float rad2 = rad*rad;
    
    //ensure that the grid contains the atom
    const AABox atombox = AABox(coords, Vector(rad+2.5));
    
    if (occ.isEmpty())
    {
        grid_info = GridInfo(atombox, grid_spacing);
        
        if (grid_info.nPoints() > max_points)
        {
            grid_info = GridInfo();
            throw SireError::unavailable_resource( QObject::tr(
                    "Unable to add the atoms as doing so would increase the number of grid points "
                    "to beyond the maximum supported size. Adding the atom at position %1, radius "
                    "%2 would create the grid %3, requiring %4 grid points. This is "
                    "greater than the maximum number of grid points supported (%5). To add "
                    "this atom, either increase the grid spacing, or increase the maximum number "
                    "of allowable grid points.")
                        .arg(coords.toString()).arg(rad)
                        .arg(grid_info.dimensions().toString())
                        .arg(grid_info.nPoints())
                        .arg(max_points), CODELOC );
        }
        
        occ = QVector<float>(grid_info.nPoints(), 0.0);
    }
    else if (not grid_info.dimensions().contains(atombox))
    {
        //see if this atom is within the mask distance (if it is masked)
        if (not mask_points.isEmpty())
        {
            bool is_masked = true;
            
            Cartesian space;
            
            for (int i=0; i<mask_points.count(); ++i)
            {
                if (space.minimumDistance(mask_points.at(i),atombox) <= mask_dist)
                {
                    is_masked = false;
                    break;
                }
            }
            
            if (is_masked)
                //this atom will be masked, so don't bother evaluating it
                return;
        }
    
        //must redimension
        Vector mincoords = grid_info.dimensions().minCoords();
        Vector maxcoords = grid_info.dimensions().maxCoords();
        
        while (mincoords.x() > atombox.minCoords().x())
        {
            mincoords.setX( mincoords.x() - grid_spacing.value() );
        }

        while (mincoords.y() > atombox.minCoords().y())
        {
            mincoords.setY( mincoords.y() - grid_spacing.value() );
        }

        while (mincoords.z() > atombox.minCoords().z())
        {
            mincoords.setZ( mincoords.z() - grid_spacing.value() );
        }
        
        while (maxcoords.x() < atombox.maxCoords().x())
        {
            maxcoords.setX( maxcoords.x() + grid_spacing.value() );
        }

        while (maxcoords.y() < atombox.maxCoords().y())
        {
            maxcoords.setY( maxcoords.y() + grid_spacing.value() );
        }

        while (maxcoords.z() < atombox.maxCoords().z())
        {
            maxcoords.setZ( maxcoords.z() + grid_spacing.value() );
        }
        
        GridInfo new_grid_info( AABox::from(mincoords,maxcoords), grid_spacing );
        
        if (new_grid_info.nPoints() > max_points)
            throw SireError::unavailable_resource( QObject::tr(
                    "Unable to add the atoms as doing so would increase the number of grid points "
                    "to beyond the maximum supported size. Adding the atom at position %1, radius "
                    "%2 would extend the grid from %3 to %4, requiring %5 grid points. This is "
                    "greater than the maximum number of grid points supported (%6). To add "
                    "this atom, either increase the grid spacing, or increase the maximum number "
                    "of allowable grid points.")
                        .arg(coords.toString()).arg(rad)
                        .arg(grid_info.dimensions().toString())
                        .arg(new_grid_info.dimensions().toString())
                        .arg(new_grid_info.nPoints())
                        .arg(max_points), CODELOC );

        occ = grid_info.redimension(occ, new_grid_info);

        grid_info = new_grid_info;
    }
    
    GridIndex idx = grid_info.pointToGridIndex(coords);
    
    const int nboxes = int(rad / grid_spacing.value()) + 2;

    const Vector center_point = grid_info.point(idx);

    QVarLengthArray<int> occupied_points;
    
    for (int i=-nboxes; i<=nboxes; ++i)
    {
        for (int j=-nboxes; j<=nboxes; ++j)
        {
            for (int k=-nboxes; k<=nboxes; ++k)
            {
                Vector grid_point = center_point + Vector(i*grid_spacing.value(),
                                                          j*grid_spacing.value(),
                                                          k*grid_spacing.value());
            
                if (Vector::distance2(grid_point,coords) <= rad2)
                {
                    int point_idx = grid_info.gridToArrayIndex(idx.i() + i,
                                                               idx.j() + j,
                                                               idx.k() + k);

                    if (point_idx < 0 or point_idx >= occ.count())
                        throw SireError::program_bug( QObject::tr(
                                "Could not add the point as dimensioning was incorrect?"),
                                    CODELOC );
                
                    occupied_points.push_back(point_idx);
                }
            }
        }
    }
    
    //we've found all of the occupied points for this atom - add them onto the grid
    float *o = occ.data();
    
    for (int i=0; i<occupied_points.count(); ++i)
    {
        if (not mask_points.isEmpty())
        {
            bool is_masked = true;
            
            for (int j=0; j<mask_points.count(); ++j)
            {
                if (Vector::distance(mask_points.at(j),
                                     grid_info.point(occupied_points.at(i))) <= mask_dist)
                {
                    is_masked = false;
                    break;
                }
            }
        
            if (not is_masked)
                o[ occupied_points.constData()[i] ] = 1;
        }
        else
            o[ occupied_points.constData()[i] ] = 1;
    }
}

void VolumeMap::evaluate(const MoleculeView &molecule, const PropertyMap &map)
{
    if (not d or molecule.isEmpty())
        return;

    //extract the coordinates and element properties from the molecule
    const AtomCoords &coords = molecule.data().property( map["coordinates"] ).asA<AtomCoords>();
    const AtomElements &elems = molecule.data().property( map["element"] ).asA<AtomElements>();
    
    if (molecule.selectedAll())
    {
        //loop over all atoms and add them to the grid
        for (int i=0; i<coords.nCutGroups(); ++i)
        {
            const Vector *ca = coords.data(CGIdx(i));
            const Element *ea = elems.data(CGIdx(i));
            
            for (int j=0; j<coords.nAtoms(CGIdx(i)); ++j)
            {
                if ((not skipLightAtoms()) or ea[j].nProtons() >= 6)
                {
                    d->addToGrid(ca[j], ea[j]);
                }
            }
        }
    }
    else
    {
        //not selected all atoms...
        const AtomSelection selected_atoms = molecule.selection();
        
        for (int i=0; i<coords.nCutGroups(); ++i)
        {
            const Vector *ca = coords.data(CGIdx(i));
            const Element *ea = elems.data(CGIdx(i));

            if (selected_atoms.selectedAll(CGIdx(i)))
            {
                for (int j=0; j<coords.nAtoms(CGIdx(i)); ++j)
                {
                    if ((not skipLightAtoms()) or ea[j].nProtons() >= 6)
                    {
                        d->addToGrid(ca[j], ea[j]);
                    }
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(CGIdx(i)))
                {
                    if ((not skipLightAtoms()) or ea[j.value()].nProtons() >= 6)
                    {
                        d->addToGrid(ca[j.value()], ea[j.value()]);
                    }
                }
            }
        }
    }
}

void VolumeMap::evaluate(const Molecules &molecules, const PropertyMap &map)
{
    if (not d)
        return;
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        evaluate(it.value(), map);
    }
}

void VolumeMap::evaluate(const GridInfo &grid, const QVector<float> &vals, qint64 n)
{
    if (not d)
        return;

    throw SireError::incomplete_code( QObject::tr(
            "Haven't yet written the code to add VolumeMaps together...!"), CODELOC );
}

void VolumeMap::cancelEvaluation()
{
    delete d;
    d = 0;
}

void VolumeMap::endEvaluation()
{
    if (not d)
        return;
    
    if (grid_info != d->grid_info)
    {
        if (occ.isEmpty())
        {
            occ = QVector<float>(d->grid_info.nPoints(), 0.0);
        }
        else
        {
            occ = grid_info.redimension(occ, d->grid_info);
        }
        
        grid_info = d->grid_info;
    }
    else if (occ.isEmpty())
    {
        occ = QVector<float>(grid_info.nPoints(), 0.0);
    }

    if (occ.count() != d->occ.count())
        throw SireError::program_bug( QObject::tr(
                "Problem finalising grid evaluation: %1 vs. %2.")
                    .arg(occ.count()).arg(d->occ.count()), CODELOC );

    switch (map_type)
    {
        case AVERAGE:
        {
            const float big_ratio = float(nsamples) / float(nsamples+1);
            const float small_ratio = float(1) / float(nsamples+1);
            
            for (int i=0; i<occ.count(); ++i)
            {
                occ[i] = big_ratio*occ[i] + small_ratio*d->occ[i];
            }
        }
        break;
        
        case MAXIMUM:
        {
            for (int i=0; i<occ.count(); ++i)
            {
                occ[i] = qMax(occ[i], d->occ[i]);
            }
        }
        break;
        
        case SUM:
        {
            for (int i=0; i<occ.count(); ++i)
            {
                occ[i] += d->occ[i];
            }
        }
        break;
    }

    nsamples += 1;

    //qint64 nsecs = d->t.nsecsElapsed();
    delete d;
    d = 0;

    //qDebug() << "VolumeMapping took" << (0.000001*nsecs) << "ms";
}

/** Add a single molecule to the map */
void VolumeMap::add(const MoleculeView &molecule, const PropertyMap &map)
{
    try
    {
        beginEvaluation();
        evaluate(molecule, map);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add two molecules to the map */
void VolumeMap::add(const MoleculeView &mol0, const MoleculeView &mol1,
                    const PropertyMap &map)
{
    try
    {
        beginEvaluation();
        evaluate(mol0, map);
        evaluate(mol1, map);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add a whole list of molecules to the map */
void VolumeMap::add(const QList<PartialMolecule> &molecules, const PropertyMap &map)
{
    try
    {
        beginEvaluation();
        foreach (const PartialMolecule &molecule, molecules)
        {
            evaluate(molecule, map);
        }
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add a set of molecules to the map */
void VolumeMap::add(const Molecules &molecules, const PropertyMap &map)
{
    if (occ.isEmpty())
        presize(molecules, map);

    try
    {
        beginEvaluation();
        evaluate(molecules, map);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add two sets of molecules to the map */
void VolumeMap::add(const Molecules &mols0, const Molecules &mols1, const PropertyMap &map)
{
    if (occ.isEmpty())
    {
        if (mols0.nMolecules() >= mols1.nMolecules())
            presize(mols0, map);
        else
            presize(mols1, map);
    }

    try
    {
        beginEvaluation();
        evaluate(mols0, map);
        evaluate(mols1, map);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add a whole list of sets of molecules to the map */
void VolumeMap::add(const QList<Molecules> &molecules, const PropertyMap &map)
{
    if (occ.isEmpty())
    {
        int biggest = -1;
        int nbiggest = 0;
        
        for (int i=0; i<molecules.count(); ++i)
        {
            if (molecules.at(i).nMolecules() > nbiggest)
            {
                biggest = i;
                nbiggest = molecules.at(i).nMolecules();
            }
        }
        
        if (biggest >= 0)
            presize( molecules.at(biggest), map );
    }

    try
    {
        beginEvaluation();
        foreach (const Molecules &mols, molecules)
        {
            evaluate(mols, map);
        }
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add a moleculegroup to the map */
void VolumeMap::add(const MoleculeGroup &molecules, const PropertyMap &map)
{
    if (occ.isEmpty())
        presize(molecules.molecules(), map);

    try
    {
        beginEvaluation();
        evaluate(molecules.molecules(), map);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/* Add two moleculegroups to the map */
void VolumeMap::add(const MoleculeGroup &mols0, const MoleculeGroup &mols1,
                    const PropertyMap &map)
{
    if (occ.isEmpty())
    {
        if (mols0.nMolecules() >= mols1.nMolecules())
            presize(mols0.molecules(), map);
        else
            presize(mols1.molecules(), map);
    }

    try
    {
        beginEvaluation();
        evaluate(mols0.molecules(), map);
        evaluate(mols1.molecules(), map);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add a whole list of molecule groups to the map */
void VolumeMap::add(const QList<MoleculeGroup> &molecules, const PropertyMap &map)
{
    if (occ.isEmpty())
    {
        int biggest = -1;
        int nbiggest = 0;
        
        for (int i=0; i<molecules.count(); ++i)
        {
            if (molecules.at(i).nMolecules() > nbiggest)
            {
                biggest = i;
                nbiggest = molecules.at(i).nMolecules();
            }
        }
        
        if (biggest >= 0)
            presize( molecules.at(biggest).molecules(), map );
    }

    try
    {
        beginEvaluation();
        foreach (const MoleculeGroup &group, molecules)
        {
            evaluate(group.molecules(),map);
        }
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add the data from the other passed volume map onto this map */
void VolumeMap::add(const VolumeMap &other)
{
    try
    {
        beginEvaluation();
        evaluate(other.gridInfo(), other.occupancy(), other.nsamples);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Add the data from the passed grid onto this map */
void VolumeMap::add(const GridInfo &gridinfo, const QVector<float> &values)
{
    try
    {
        beginEvaluation();
        evaluate(gridinfo, values, this->nsamples);
        endEvaluation();
    }
    catch(...)
    {
        cancelEvaluation();
        throw;
    }
}

/** Return whether or not this volume map is masked */
bool VolumeMap::isMasked() const
{
    return not mask_points.isEmpty();
}

/** Return the mask distance. Grid points are only evaluated if they are
    within this distance of any of the masking points */
Length VolumeMap::maskDistance() const
{
    return Length(mask_dist);
}

/** Return all of the masking points. Grid points are only evaluated if
    they are within the mask distance of any of these points */
QVector<Vector> VolumeMap::maskPoints() const
{
    return mask_points;
}

/** Internal function used to clear points that lie outside the current mask */
void VolumeMap::applyMask()
{
    if (mask_points.isEmpty())
        return;

    for (int i=0; i<grid_info.nPoints(); ++i)
    {
        Vector c = grid_info.point(i);
        
        bool should_be_masked = true;
        
        for (int j=0; j<mask_points.count(); ++j)
        {
            Vector m = mask_points.at(j);
            
            if (Vector::distance(c,m) <= mask_dist)
            {
                should_be_masked = false;
                break;
            }
        }
        
        if (should_be_masked)
            occ[i] = 0;
    }
}

/** Set the mask such that grid points are only evaluated if they are within
    'dist' distance of point 'point'. If 'clear_points' is true (default), then this
    will clear any points that are outside this mask */
void VolumeMap::setMaskWithinDistance(Length dist, const Vector &point, bool clear_points)
{
    mask_points.clear();
    mask_points.append(point);
    mask_dist = dist.value();
    
    if (mask_dist < 0)
        mask_dist = 0;
    
    if (clear_points)
        applyMask();
}

/** Set the mask such that grid points are only evaluated if they are within
    distance 'dist' of any atom in the view 'molecule'. If 'clear_points' is true
    (default) then this will clear any points that are outside the mask */
void VolumeMap::setMaskWithinDistance(Length dist, const MoleculeView &molecule,
                                      bool clear_points, const PropertyMap &map)
{
    const AtomCoords &coords = molecule.data().property( map["coordinates"] )
                                              .asA<AtomCoords>();
    
    mask_points.clear();
    mask_dist = dist.value();

    if (mask_dist < 0)
        mask_dist = 0;
    
    if (molecule.selectedAll())
    {
        for (int i=0; i<coords.nCutGroups(); ++i)
        {
            const Vector *c = coords.data(CGIdx(i));
        
            for (int j=0; j<coords.nAtoms(CGIdx(i)); ++j)
            {
                mask_points.append( c[j] );
            }
        }
    }
    else
    {
        AtomSelection selected_atoms = molecule.selection();
        
        for (int i=0; i<coords.nCutGroups(); ++i)
        {
            CGIdx ci(i);
        
            if (selected_atoms.selected(ci))
            {
                const Vector *c = coords.data(ci);
            
                if (selected_atoms.selectedAll(ci))
                {
                    for (int j=0; j<coords.nAtoms(ci); ++j)
                    {
                        mask_points.append( c[j] );
                    }
                }
                else
                {
                    foreach (Index j, selected_atoms.selectedAtoms(ci))
                    {
                        mask_points.append( c[j] );
                    }
                }
            }
        }
    }
    
    if (clear_points)
        applyMask();
}

/** Set the mask such that grid points are only evaluated if they are within
    distance 'dist' of any atom in the view 'molecule'. If 'clear_points' is true
    (default) then this will clear any points that are outside the mask */
void VolumeMap::setMaskWithinDistance(Length dist, const MoleculeView &molecule,
                                      const PropertyMap &map)
{
    setMaskWithinDistance(dist, molecule, true, map);
}

/** Clear the set of mask points and mask distance */
void VolumeMap::clearMask()
{
    mask_points.clear();
    mask_dist = 0;
}

