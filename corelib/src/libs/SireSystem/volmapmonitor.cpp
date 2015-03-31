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

#include "volmapmonitor.h"

#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomcoords.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/mgnum.h"

#include "SireVol/gridinfo.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QElapsedTimer>

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireVol;
using namespace SireStream;
using namespace SireUnits::Dimension;
using namespace SireUnits;

static const RegisterMetaType<VolMapMonitor> r_mon;

QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, const VolMapMonitor &mon)
{
    writeHeader(ds, r_mon, 1);

    SharedDataStream sds(ds);
    
    sds << mon.molgroup << mon.propmap << mon.volmap << static_cast<const SystemMonitor&>(mon);
    
    return ds;
}

QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, VolMapMonitor &mon)
{
    VersionID v = readHeader(ds, r_mon);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> mon.molgroup >> mon.propmap >> mon.volmap >> static_cast<SystemMonitor&>(mon);
    }
    else
        throw version_error(v, "1", r_mon, CODELOC);
    
    return ds;
}

static const VolumeMap::MapType map_type = VolumeMap::AVERAGE;
static const VolumeMap::FillType fill_type = VolumeMap::VDW_RADIUS;

/** Null constructor */
VolMapMonitor::VolMapMonitor()
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                volmap(map_type, fill_type)
{}

/** Construct specifying the grid spacing */
VolMapMonitor::VolMapMonitor(const SireUnits::Dimension::Length &grid_spacing)
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                volmap(grid_spacing, map_type, fill_type)
{}

/** Construct, specifying the molecule group to be monitored */
VolMapMonitor::VolMapMonitor(const MoleculeGroup &group, const PropertyMap &map)
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                molgroup(group), propmap(map), volmap(map_type, fill_type)
{}

/** Construct, specifying the molecule group to be monitored and the grid spacing 
    for the occupancy grid */
VolMapMonitor::VolMapMonitor(const MoleculeGroup &group,
                             const SireUnits::Dimension::Length &grid_spacing,
                             const PropertyMap &map)
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                molgroup(group), propmap(map), volmap(grid_spacing,fill_type,map_type)
{}

/** Construct, specifying the molecule group to be monitored
    and whether or not to ignore light atoms */
VolMapMonitor::VolMapMonitor(const MoleculeGroup &group,
                             bool skip_light_atoms,
                             const PropertyMap &map)
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                molgroup(group), propmap(map), volmap(fill_type,map_type)
{
    volmap.setSkipLightAtoms(skip_light_atoms);
}

/** Construct, specifying the molecule group to be monitored, the grid spacing
    for the occupancy grid, and whether or not to ignore light atoms */
VolMapMonitor::VolMapMonitor(const MoleculeGroup &group,
                             const SireUnits::Dimension::Length &grid_spacing,
                             bool skip_light_atoms,
                             const PropertyMap &map)
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                molgroup(group), propmap(map), volmap(grid_spacing,map_type,fill_type)
{
    volmap.setSkipLightAtoms(skip_light_atoms);
}

/** Copy constructor */
VolMapMonitor::VolMapMonitor(const VolMapMonitor &other)
              : ConcreteProperty<VolMapMonitor,SystemMonitor>(),
                molgroup(other.molgroup), propmap(other.propmap), volmap(other.volmap)
{}

/** Destructor */
VolMapMonitor::~VolMapMonitor()
{}

/** Copy assignment operator */
VolMapMonitor& VolMapMonitor::operator=(const VolMapMonitor &other)
{
    if (this != &other)
    {
        molgroup = other.molgroup;
        volmap = other.volmap;
        propmap = other.propmap;
    }
    
    return *this;
}

/** Comparison operator */
bool VolMapMonitor::operator==(const VolMapMonitor &other) const
{
    return molgroup == other.molgroup and
           propmap == other.propmap and
           volmap == other.volmap;
}

/** Comparison operator */
bool VolMapMonitor::operator!=(const VolMapMonitor &other) const
{
    return not operator==(other);
}

const char* VolMapMonitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VolMapMonitor>() );
}

const char* VolMapMonitor::what() const
{
    return VolMapMonitor::typeName();
}

QString VolMapMonitor::toString() const
{
    return QObject::tr("VolMapMonitor( group() == %1, gridInfo() == %2, nSamples() == %3 )")
                .arg(group().toString())
                .arg(gridInfo().toString())
                .arg(nSamples());
}

/** Return the molecule group whose atoms are being monitored */
const MoleculeGroup& VolMapMonitor::group() const
{
    return molgroup.read();
}

/** Set the molecule group being monitored (together with the property map 
    used to get the required properties from that group) */
void VolMapMonitor::setGroup(const MoleculeGroup &group, const PropertyMap &map)
{
    clearStatistics();
    molgroup = group;
    propmap = map;
}

/** Return whether or not we are skipping light atoms (e.g. hydrogen) */
bool VolMapMonitor::skippingLightAtoms() const
{
    return volmap.skipLightAtoms();
}

/** Turn on skipping light atoms */
void VolMapMonitor::setSkipLightAtoms(bool on)
{
    volmap.setSkipLightAtoms(on);
}

/** Return the spacing of the grid */
Length VolMapMonitor::gridSpacing() const
{
    return volmap.gridSpacing();
}

/** Set the grid spacing on which we monitor the space */
void VolMapMonitor::setGridSpacing(const Length &grid_spacing)
{
    volmap.setGridSpacing(grid_spacing);
}

/** Return the grid dimensions */
GridInfo VolMapMonitor::gridInfo() const
{
    return volmap.gridInfo();
}

/** Return the property map used to find the properties needed by this monitor */
PropertyMap VolMapMonitor::map() const
{
    return propmap;
}

/** Set the property map to be used by this monitor */
void VolMapMonitor::setPropertyMap(const PropertyMap &m)
{
    propmap = m;
}

/** Return the average occupancy. This is a linear array that can be
    accessed using the accompanying GridInfo returned by gridInfo() */
QVector<float> VolMapMonitor::averageOccupancy() const
{
    return volmap.occupancy();
}

/** Return the number of samples recorded so far */
qint64 VolMapMonitor::nSamples() const
{
    return volmap.nSamples();
}

/** Clear all statistics */
void VolMapMonitor::clearStatistics()
{
    volmap.clear();
}

/** Return the actual volume map */
VolumeMap VolMapMonitor::volumeMap() const
{
    return volmap;
}

/** Monitor the system */
void VolMapMonitor::monitor(System &system)
{
    MolGroupPtr group_to_monitor;
    
    //now update the molecule group, if needed
    if (system.contains(group().number()))
    {
        group_to_monitor = system[group().number()];
    }
    else
    {
        group_to_monitor = molgroup;
        group_to_monitor.edit().update(system.molecules());
    }
    
    if (group_to_monitor.read().isEmpty())
        return;

    volmap.add(group_to_monitor.read(), propmap);
}
